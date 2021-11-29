#!/bin/bash

rsync -avzL --files-from=Spectra_to_Download_results.csv_1.txt rsync://data.sdss.org/dr16 dr16