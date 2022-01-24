from cmath import log
import json
import logging
import coloredlogs
from dataclasses import asdict
from flask import Flask, render_template, request, redirect, url_for
from rratsolve import rratsolve


app = Flask(__name__)
last_result = None


def set_last_result(r):
    global last_result
    last_result = r


def parse_toas(toas_str):
    seq = toas_str.replace('\r', '').replace(',', '').split('\n')
    seq = filter(len, seq)
    return list(map(float, seq))


@app.route('/result.json')
def result_json():
    return json.dumps(asdict(last_result), indent=4)


@app.route('/solve')
def solve():
    toas = request.args.get('toas')
    toa_uncertainty = request.args.get('uncertainty')
    T = parse_toas(toas)
    u = float(toa_uncertainty)
    set_last_result(rratsolve(T, u, max_grid_size=30_000_000))
    return render_template('solve.html')


@app.route('/')
def mainpage():
    return render_template('index.html')


def main():
    app.run(debug=True)


if __name__ == '__main__':
    coloredlogs.install(logger=logging.getLogger('rratsolve'), level='DEBUG')
    main()