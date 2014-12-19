#!/bin/bash
ls -a -1 *.sam > samfiles
sed -i 's/.sam//g' samfiles
