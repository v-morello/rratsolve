from cmath import log
import json
import logging
import coloredlogs
from flask import Flask, render_template, request, redirect, url_for
from rratsolve import rratsolve
from rratsolve.core import format_uncertain_quantity

app = Flask(__name__)


def parse_toas(toas_str):
    print(f"toas_str = {toas_str!r}")
    seq = toas_str.replace('\r', '').replace(',', '').split('\n')
    seq = filter(len, seq)
    return list(map(float, seq))


@app.route('/solve')
def call_solver():
    toas = request.args.get('toas')
    toa_uncertainty = request.args.get('uncertainty')
    T = parse_toas(toas)
    u = float(toa_uncertainty)
    period, perr = rratsolve(T, u)
    period_str = format_uncertain_quantity(period, perr)

    toas_formatted = "\n    " + "\n    ".join(map(str, T))
    print(toas_formatted)
    return render_template('solve.html', period=period_str, toas=f"{toas_formatted}", toa_uncertainty=u)


@app.route('/')
def mainpage():
    return render_template('index.html')


def main():
    app.run(debug=True)


if __name__ == '__main__':
    coloredlogs.install(logger=logging.getLogger('rratsolve'), level='DEBUG')
    main()