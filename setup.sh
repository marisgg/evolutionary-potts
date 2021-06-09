#!/bin/bash
sudo apt update
sudo apt upgrade -y
sudo apt install nodejs npm -y
npm install
sudo apt install python3 python3-pip -y
pip install -r requirements.txt