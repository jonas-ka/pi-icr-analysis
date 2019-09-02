# -*- coding: utf-8 -*-
"""
Created on 17 July 2019

@author: Jonas Karthein
@contact: jonas.karthein@cern.ch
@license: MIT license
"""

import sys, os, glob, json, csv
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize



class Freq_ratio():
    # class to calculate the frequency ratio for a mass measurement (MM), optimized for ISOLTRAP PI-ICR data
    def __init__(self):
        self.degree = 3
        self.param_start = [1]*(self.degree+2)
        self.param_result = []
        self.y1df = {}
        self.y1names = []
        self.y1 = []
        self.y1err = []
        self.y2df = {}
        self.y2names = []
        self.y2 = []
        self.y2err = []
        self.x1df = {}
        self.x1helpdf = {}
        self.x1names = []
        self.x1 = []
        self.x2df = {}
        self.x2helpdf = {}
        self.x2names = []
        self.x2 = []
        self.x2folder = ''      # time stamps isotope 2
        self.isotopes = []
        self.dirname = ''
        self.run_folder = ''
        self.ratio = 0
        self.ratio_unc = 0
        self.red_chi_sq = 1
        self.covar = []
        self.correlation = 0
        self.ratio_dict = {}
        self.weighted_avg = 0
        self.total_unc = 0
        self.birge_ratio = 0
        self.merged_dict = {}
        self.lin_extrapolation = False
        self.mode = ''

        self.test = False

    def ratio_sim_fit(self, isotopes, cyc_freq_1_time_delta, cyc_freq_1, cyc_freq_1_unc, cyc_freq_2_time_delta, cyc_freq_2, cyc_freq_2_unc):
        """"
        Main function perform simultaneous polynomial fits for all for the data set possible polynomial degrees. Calculates red. chi square for each fit and returns the fit results for the one with a red. chi square clostest to 1.

        :param isotopes: list of str of isotopes
        :param cyc_freq_1_time_delta: list of time deltas in minutes to first time for isotope 1
        :param cyc_freq_1: list of cyclotron frequencies for isotope 1
        :param cyc_freq_1_unc: list of cyclotron frequency uncertainties for isotope 1
        :param cyc_freq_2_time_delta: list of time deltas in minutes to first time for isotope 2
        :param cyc_freq_2: list of cyclotron frequencies for isotope 2
        :param cyc_freq_2_unc: list of cyclotron frequency uncertainties for isotope 2

        :return param1: dict containing the fit information for isotope 1
        :return param2: dict containing the fit information for isotope 2
        :return ratio: frequency ratio deriving from best fit
        :return ratio_unc: uncertainty of frequency ratio deriving from best fit
        :return red_chi_sq: reduced chi square for best fit
        """
        self.x1 = cyc_freq_1_time_delta
        self.x2 = cyc_freq_2_time_delta
        self.y1 = cyc_freq_1
        self.y1err = cyc_freq_1_unc
        self.y2 = cyc_freq_2
        self.y2err = cyc_freq_2_unc
        self.isotopes = isotopes
        temp_result_dict = {}
        temp_result_dict2 = {}

        for deg in range(2,min(len(self.y1),len(self.y2))+2):     # calculates ratio for many polynomial degrees, as much as data allows
            try:
                self.degree = deg
                self.get_start_param()
                if len(self.y2) > 1:    # make sure there are at least two points to fit
                    print('Poly-degree: ', self.degree)
                    self.fit()
                    self.red_chi_sq_calc()
                    print('\nRed.Chi.Sq.: ', self.red_chi_sq)
                    self.covar = self.covar * self.red_chi_sq
                    self.param_result_err = np.array([np.sqrt(np.absolute(self.covar[i, i])) for i in range(len(self.param_result))])
                    self.correlation = self.covar[len(self.param_result)-2, len(self.param_result)-1]       # calculates the corellation
                    print('Corellation: ', self.correlation)
                else:       # if less than two points in data set, then a linear extrapolation is performed
                    self.lin_extrapolation = True
                    self.fit()
                    self.param_result_err = np.sqrt(np.diag(self.covar))
                # calculates the frequency ratio from the fit information
                self.ratio_calc()
                self.lin_extrapolation = False

                temp_result_dict['{}'.format(self.degree)] = {'param1': [float(i) for i in list(self.param_result[:-1])],
                                                              'param2': [float(i)/float(self.param_result[-1]) for i in list(self.param_result[:-1])],
                                                              'ratio': self.ratio,
                                                              'ratio_unc': self.ratio_unc,
                                                              'red_chi_sq': self.red_chi_sq}
                temp_result_dict2['{}'.format(self.degree)] = abs(self.red_chi_sq-1.0)      # allowing to find the red_chi_sq closest to 1 by a minimizing function
            except:     # if one fit fails, program stops
                break

        try:
            best_chi_sq = min(temp_result_dict2, key=temp_result_dict2.get)
        except:
            best_chi_sq = -1

        return(temp_result_dict[best_chi_sq]['param1'], temp_result_dict[best_chi_sq]['param2'], temp_result_dict[best_chi_sq]['ratio'], temp_result_dict[best_chi_sq]['ratio_unc'], temp_result_dict[best_chi_sq]['red_chi_sq'])

    def get_start_param(self):
        """Sets starting parameter. If fitting fails, try to adjust the starting parameter"""
        self.param_start = [1]*(self.degree+1)
        self.param_start[-2] = self.y1[0]   # linear offset of fit
        self.param_start[-3] = 0.0001           # helps p=2 fit to converge
        #     self.param_start[-3] = 0.1           # helps p=2 fit to converge
        #     print(self.y2, self.y1)
        print(self.param_start)

    def get_data(self):
        """Get frequency and timestamps"""
        # ___________________________________________
        # get frequency data for isotope 1
        self.y1df = pd.read_csv(self.y1folder)
        self.y1names = self.y1df.iloc[:, 0].tolist()
        self.y1df = self.y1df.sort_values('File')
        self.y1df = self.y1df.reset_index(drop=True)
        self.y1 = self.y1df.iloc[:, 1].tolist()
        self.y1err = self.y1df.iloc[:, 2].tolist()

        # ___________________________________________
        # get frequency data for isotope 2
        self.y2df = pd.read_csv(self.y2folder)
        self.y2df = self.y2df.sort_values('File')
        self.y2df = self.y2df.reset_index(drop=True)
        self.y2names = self.y2df.iloc[:, 0].tolist()
        self.y2 = self.y2df.iloc[:, 1].tolist()
        self.y2err = self.y2df.iloc[:, 2].tolist()
        # ___________________________________________
        # get time stamps for each file for isotope 1
        self.x1df = pd.read_csv(self.x1folder,
                                names=['names', 'date_begin', 'time_begin', 'date_end', 'time_end'],
                                parse_dates = [[1,2], [3,4]])
        self.x1df = self.x1df.sort_values('names')
        self.x1df = self.x1df.reset_index(drop=True)
        self.x1df['timedelta'] = self.x1df['date_end_time_end'] - self.x1df['date_begin_time_begin']
        self.x1names = self.x1df['names'].tolist()
        # calculate the median time stamp for each file also considering the count rate
        self.x1df['median_counts'] = self.get_timestamp_position(self.x1helpfolder, self.x1names)
        self.x1df['median_time'] = (self.x1df['date_begin_time_begin'] +
                                    self.x1df['timedelta'] * self.x1df['median_counts'])


        # ___________________________________________
        # get time stamps for each file for isotope 2
        self.x2df = pd.read_csv(self.x2folder,
                                names=['names', 'date_begin', 'time_begin', 'date_end', 'time_end'],
                                parse_dates = [[1,2], [3,4]])
        self.x2df = self.x2df.sort_values('names')
        self.x2df = self.x2df.reset_index(drop=True)

        self.x2df['timedelta'] = self.x2df['date_end_time_end'] - self.x2df['date_begin_time_begin']
        self.x2names = self.x2df['names'].tolist()
        # calculate the median time stamp for each file also considering the count rate
        self.x2df['median_counts'] = self.get_timestamp_position(self.x2helpfolder, self.x2names)
        self.x2df['median_time'] = (self.x2df['date_begin_time_begin'] +
                                    self.x2df['timedelta'] * self.x2df['median_counts'])
        # the x values for the frequency information is converted into minutes after(/before) the first
        # x1 frequency point started. This had to be changed since the fit wasn't very stable for x values
        # having the correct time converted to seconds since the values were very large (e.g. 1498777404).
        # The following statement reduces the x value to minutes in the 2-3 digit range. Important here is,
        # that x1 and x2 are normalized to the same starting time of MIN(x1[0], x2[0])
        self.x1 = [round(i.total_seconds()/60, 2) for i in (self.x1df['median_time'] - min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).tolist()]
        self.x2 = [round(i.total_seconds()/60, 2) for i in (self.x2df['median_time'] - min([self.x1df['median_time'][0], self.x2df['median_time'][0]])).tolist()]


    def get_timestamp_position(self, folder, x1_or_x2_names):
        """calculates the mean position of ions to arrive during the measurement time --> if beam gate was changed during the run, the average timestamp is not in the middle of the measurement"""
        median_counts = []
        for i in x1_or_x2_names:
            tmp = folder+i+'_p1_spot_positions.csv'
            self.x1helpdf = pd.read_csv(tmp, names=['event', 'x', 'y', 'time_stamp_us'])
            # we chose the median over the mean to be less sensitive to outliers in case of a drop in count rate
            median_counts.append(float(self.x1helpdf.iloc[:, 0].median() / self.x1helpdf.iloc[:, 0].max()))
        return(median_counts)


    def fit(self):
        """Actual simultaneous polynomial fitting routine. Executed five times for stability reasons."""
        if self.lin_extrapolation == False:
            # fit ; timestampt value is reduced since resulting numbers were too large
            iter1, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, self.param_start,
                                                            args=(np.array(self.x1),
                                                                  np.array(self.x2),
                                                                  np.array(self.y1),
                                                                  np.array(self.y2),
                                                                  np.array(self.y1err),
                                                                  np.array(self.y2err)), full_output=True)#, maxfev=1000000
            iter2, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter1,
                                                            args=(np.array(self.x1),
                                                                  np.array(self.x2),
                                                                  np.array(self.y1),
                                                                  np.array(self.y2),
                                                                  np.array(self.y1err),
                                                                  np.array(self.y2err)), full_output=True)
            iter3, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter2,
                                                            args=(np.array(self.x1),
                                                                  np.array(self.x2),
                                                                  np.array(self.y1),
                                                                  np.array(self.y2),
                                                                  np.array(self.y1err),
                                                                  np.array(self.y2err)), full_output=True)
            iter4, covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter3,
                                                            args=(np.array(self.x1),
                                                                  np.array(self.x2),
                                                                  np.array(self.y1),
                                                                  np.array(self.y2),
                                                                  np.array(self.y1err),
                                                                  np.array(self.y2err)), full_output=True)
            self.param_result, self.covar, infodict, mesg, ier = scipy.optimize.leastsq(self.simultaneous_residual, iter4,
                                                            args=(np.array(self.x1),
                                                                  np.array(self.x2),
                                                                  np.array(self.y1),
                                                                  np.array(self.y2),
                                                                  np.array(self.y1err),
                                                                  np.array(self.y2err)), full_output=True)
            if ier not in [1,2,3,4]:
                sys.exit('Fit did not converge')

            # print self.param_result
            # print self.covar

        else:       # linear extrapolation
            init_vals = [((max(self.y1)-min(self.y1))
                          /(max(self.x1)-min(self.x1))),self.y1[0]]
            fit_1, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=init_vals)
            fit_2, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_1)
            fit_3, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_2)
            fit_4, pcov = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_3)
            self.param_result, self.covar = scipy.optimize.curve_fit(self.linear,
                                       self.x1, self.y1, sigma=self.y1err,
                                       absolute_sigma=True,
                                       p0=fit_4)
            print(self.param_result)

    def fit_function(self, x, param):
        """
        Returns a polynomial function value of degree len(param)-1 evaluated at x

        :param x: time delta value
        :param param: list of polynomial parameters

        :return: polynomial function value of degree len(param)-1 evaluated at x
        """
        return(np.polyval(list(param), x))

    def linear(self, x, m, n):
        """
        Returns a linear function value evaluated at x

        :param x: time delta value
        :param m: slope
        :param n: frequency axis intercept

        :return: linear function value evaluated at x
        """
        return(m*x+n)

    def residual(self, p, x, y, yerr):
        """
        Returns residual, which square is minimized later.

        :param p: list of polynomial parameters
        :param x: list of time delta values
        :param y: list of cyclotron frequency values
        :param yerr: list of cyclotron frequency uncertainty values

        :return: residual, which square is minimized later
        """
        return((self.fit_function(x, p) - y )/yerr)

    def simultaneous_residual(self, param_start, x1, x2, y1, y2, y1err, y2err):
        """
        Returns residual for both data sets.

        :param param_start: list of polynomial start parameters for fit
        :param x1: list of time delta values of isotope 1
        :param y1: list of cyclotron frequency values of isotope 1
        :param y1err: list of cyclotron frequency uncertainty values of isotope 1
        :param x2: list of time delta values of isotope 2
        :param y2: list of cyclotron frequency values of isotope 2
        :param y2err: list of cyclotron frequency uncertainty values of isotope 2

        :return: total residual
        """

        p_s = param_start
        p1 = p_s[:-1]  # Poly1 = R*Poly2
        p2 = [float(i)/p_s[-1] for i in p_s[:-1]]  # Poly1 = R*Poly2
        res1 = self.residual(p1, x1, y1, y1err)
        res2 = self.residual(p2, x2, y2, y2err)
        return(np.concatenate((res1, res2)))

    def final_residual_1(self, param_start, x1, y1, y1err):
        """
        Returns final residual for isotope 1.

        :param param_start: list of final polynomial parameters of fit
        :param x1: list of time delta values of isotope 1
        :param y1: list of cyclotron frequency values of isotope 1
        :param y1err: list of cyclotron frequency uncertainty values of isotope 1

        :return: final residual for isotope 1
        """
        if self.test:
            p_s = [param_start['param{}'.format(i)] for i in [len(param_start) - 1 - j for j in range(len(param_start))]]
        else:
            p_s = param_start
        p1 = p_s[:-1]  # Poly1 = R*Poly2
        res1 = self.residual(p1, x1, y1, y1err)
        return(res1)

    def final_residual_2(self, param_start, x2, y2, y2err):
        """
        Returns final residual for isotope 2.

        :param param_start: list of final polynomial parameters of fit
        :param x2: list of time delta values of isotope 2
        :param y2: list of cyclotron frequency values of isotope 2
        :param y2err: list of cyclotron frequency uncertainty values of isotope 2

        :return: final residual for isotope 2
        """
        if self.test:
            p_s = [param_start['param{}'.format(i)] for i in [len(param_start) - 1 - j for j in range(len(param_start))]]
        else:
            p_s = param_start
        p2 = [float(i)/p_s[-1] for i in p_s[:-1]]  # Poly1 = R*Poly2
        res2 = self.residual(p2, x2, y2, y2err)
        return(res2)

    def red_chi_sq_calc(self):
        """Calculates reduced chi square for fit."""
        sum_res_sq_1 = 0
        sum_res_sq_2 = 0
        for i in range(len(self.x1)):
            sum_res_sq_1 += (self.final_residual_1(self.param_result, self.x1[i], self.y1[i], self.y1err[i]))**2
        for i in range(len(self.x2)):
            sum_res_sq_2 += (self.final_residual_2(self.param_result, self.x2[i], self.y2[i], self.y2err[i]))**2
        self.red_chi_sq = (sum_res_sq_1+sum_res_sq_2)/(len(self.x1)+len(self.x2)-self.degree-1)

    def ratio_calc(self):
        """Calculates the ratio in the middle of the range of interest."""
        if self.lin_extrapolation == False:
            self.ratio = self.param_result[-1]
            self.ratio_unc = self.param_result_err[-1]
        else:
            fit_at_y2 = self.linear(self.x2[0], self.param_result[0], self.param_result[1])
            self.ratio = fit_at_y2 / self.y2[0]

            unc_lin = np.sqrt( ( self.param_result_err[0] * self.x2[0] )**2
                               + self.param_result_err[1] ** 2 )
            self.ratio_unc = np.sqrt( ( unc_lin / self.y2[0] )**2
                                      + ( -( fit_at_y2 * self.y2err[0] )
                                          / (self.y2[0])**2 )**2 )
            self.degree = 'linear'
            print('Degree:              ', self.degree)

        print('Ratio fit parameter: ', self.ratio, '+/-', self.ratio_unc)
