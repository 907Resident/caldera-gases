# Minimal Start Up - Web App
from flask import Flask
app = Flask(__name__)

@app.route("/")
def hello():
    return "Hello World!\n"