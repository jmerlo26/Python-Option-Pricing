from urllib.request import urlopen as uReq
from bs4 import BeautifulSoup as soup 
import html2text
import re


import numpy as np
from scipy.stats import norm

from datetime import datetime
from dateutil.parser import parse

import pprint
pp = pprint.PrettyPrinter(indent=4)

from yahoo_fin.stock_info import *
from yahoo_fin import options


def getHTMLSoup(url):
	uClient = uReq(url)
	page_html = uClient.read()
	uClient.close()
	return soup(page_html, "html.parser")

def getInterestRate():
	url = 'https://www.federalreserve.gov/releases/h15/' 
	soup = getHTMLSoup(url) #get HTML from url
	rateHTML = soup.body.find("div", id = "content").table #Parse the HTML for the rate
	tableRows = rateHTML.find_all('tr')
	for row in tableRows:
		if 'Federal funds (effective)' in row.getText(): #find key word
			return float(row.getText().split()[-1])
	return None # Returns None if an error occured or the rate was unable to be found


def getOptions(ticker):
	dates = options.get_expiration_dates(ticker) # get available expiration dates
	calls = {}
	puts = {}
	for date in dates: # pull and format data
		if (parse(date) - datetime.now()).days > 60 and (parse(date) - datetime.now()).days < 150: 
			print(date)
			calls[date] = options.get_options_chain(ticker)['calls'].set_index('Strike')
			puts[date] = options.get_options_chain(ticker)['puts'].set_index('Strike')

	return {'calls': calls, 'puts' : puts}
		
def getPrice(ticker):
	lastDayData = get_data(ticker, start_date = None, end_date = None, index_as_date = True).tail(1) # gets the most recent price 
	return lastDayData['close'].values[0]

def getPriceSeries(ticker):
	return get_data(ticker, start_date = None, end_date = None, index_as_date = True)['close'].values # gets all prices as a series

def N(z):
	return norm.cdf(z,0.0,1.0)

def NPrime(z):
	return (1/np.sqrt(2*np.pi) )* np.exp(-(z**2) / 2)

def D1(stockPrice, strike, div, rate,  sigma, T):
	return ( np.log(stockPrice / strike) + (rate - div + (sigma**2 / 2) ) * T ) / (sigma * np.sqrt(T))

def D2(stockPrice, strike, div, rate,  sigma, T):
	return (np.log(stockPrice / strike) + (rate - div - ( sigma ** 2)/2 ) * T) / (sigma * np.sqrt(T))

def BlackScholesCall(stockPrice, strike, div, rate, sigma, T):
	d1 = D1(stockPrice, strike, div,rate,sigma, T)
	d2 = D2(stockPrice, strike, div,rate, sigma, T)

	return stockPrice *np.exp(-div * T)* N(d1) - strike * np.exp(-rate * T) * N(d2)

def BlackScholesPut(stockPrice, strike, div, rate, sigma,T):
	D = np.exp(-rate * T) #Discount Factor
	F = stockPrice / D #Forward Price
	return BlackScholesCall(stockPrice, strike, div, rate, sigma, T) - (D * (F-strike)) # Put call parity



def Vega(stockPrice, strike, div, rate, sigma, T):
	return stockPrice * NPrime(D1(stockPrice, strike, div, rate, sigma, T)) * np.sqrt(T)

def IVCall(sigma, optionPrice, stockPrice, strike, div, rate, T,i):
	if i > 1000: # Function not converging to a root
		raise Exception('Not Converging')
	estimate = BlackScholesCall(stockPrice, strike, div, rate, sigma,T ) #curent price with current iv guess
	vega = Vega(stockPrice, strike, div, rate, sigma, T)
	if vega == 0:#prevent div 0 error
		raise Exception('div 0')

	newSigma = sigma + ( ( optionPrice - estimate) / vega ) #adjust iv 

	if estimate - optionPrice < 0.001 and estimate - optionPrice > -0.001: # current price if within tolerance
		return newSigma
	else:
		return IVCall(newSigma, optionPrice, stockPrice, strike, div,rate, T, i+1) # recursively solve

def IVPut(sigma, optionPrice, stockPrice, strike,div, rate, T, i):
	if i > 1000: # Function not converging to a root
		raise Exception('Not Converging')
	estimate = BlackScholesPut(stockPrice, strike, div, rate, sigma,T ) #curent price with current iv guess
	vega = Vega(stockPrice, strike,div, rate, sigma, T)
	if vega == 0: #prevent div 0 error
		raise Exception('div 0')
	newSigma = sigma - ( (estimate - optionPrice) / vega ) #adjust iv

	if estimate - optionPrice < 0.001 and estimate - optionPrice > -0.001: # return if within tolerance
		return newSigma
	else:
		return IVPut(newSigma, optionPrice, stockPrice, strike,div, rate, T, i+1) #recursively solve

def BlackScholesIV(optionPrice, stockPrice, strike, div, rate, T, Call = True):
	sigma = 0.5 # initial guess
	if Call:
		return IVCall(sigma, optionPrice, stockPrice, strike, div, rate, T,0) #solve fo call
	else: 
		return IVPut(sigma, optionPrice, stockPrice, strike, div, rate, T,0)#solve for put


		
