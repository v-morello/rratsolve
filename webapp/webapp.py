import json
import logging
import uuid
import coloredlogs
from dataclasses import asdict
from flask import Flask, render_template, request, redirect, url_for
from rratsolve import rratsolve


app = Flask(__name__)

# Dictionary {unique_id (str): result (rratsolve.Result)}
results_cache = {}


def parse_toas(toas_str):
    seq = toas_str.replace('\r', '').replace(',', '').split('\n')
    seq = filter(len, seq)
    return list(map(float, seq))


@app.route('/result/<result_id>')
def result_json(result_id=""):
    output = json.dumps(asdict(results_cache[result_id]), indent=4)
    # We can now delete the Result object from the cache
    del results_cache[result_id]
    return output


@app.route('/results')
def solve():
    toas = request.args.get('toas')
    toa_uncertainty = request.args.get('uncertainty')
    T = parse_toas(toas)
    u = float(toa_uncertainty)

    # Store result until it is fetched by the results webpage
    result_id = str(uuid.uuid4())
    results_cache[result_id] = rratsolve(T, u, max_grid_size=30_000_000)
    return render_template('results.html', result_id=result_id)


@app.route('/')
def mainpage():
    return render_template('index.html')


def main():
    app.run(debug=True)


if __name__ == '__main__':
    coloredlogs.install(logger=logging.getLogger('rratsolve'), level='DEBUG')
    main()