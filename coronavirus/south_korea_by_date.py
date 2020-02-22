'''
Created on Feb 22, 2020

@author: ywkim
'''
import datetime
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from scipy.optimize import curve_fit


def func(x, a, b, c):
    return a * np.exp(b * x) + c

if __name__ == '__main__':
    input_file = 'sk_positive_by_date.tsv'
    df = pd.read_csv(input_file, sep='\t')
    df['Date'] = pd.to_datetime(df['Date'])
    df.columns = ['Date', 'Confirmed Cases']

    padding = datetime.timedelta(days=1)
    max_padding = datetime.timedelta(days=10)
    fig, ax = plt.subplots()
    ax = sns.scatterplot(df['Date'], df['Confirmed Cases'])

    x = mdates.date2num(df['Date'].values) - 737443
    y = df['Confirmed Cases'].to_list()

    # yn = y + 0.2 * np.random.normal(size=len(x))
    # popt, pcov = curve_fit(func, x, yn)
    popt, pcov = curve_fit(func, x, y)

    fig.autofmt_xdate()

    ax.plot(x+737443, func(x, *popt), 'r-')
    ax.set_xlim([df['Date'].to_list()[0] - padding, df['Date'].to_list()[-1] + max_padding])
    plt.savefig('fitted_curve.png', dpi=300)

    x_future = np.asfarray(range(35, 41, 1))
    y_pred = func(x_future, *popt)
    ax.plot(x_future + 737443, y_pred, 'bo')

    for i, txt in enumerate(list(y_pred)):
        ax.annotate(int(txt), (x_future[i]+737444, y_pred[i]-500), size=8)

    plt.savefig('fitted_curve_with_prediction.png', dpi=300)
