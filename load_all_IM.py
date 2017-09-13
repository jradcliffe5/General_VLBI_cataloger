#!/usr/bin/env ParselTongue
import AIPS
from Wizardry.AIPSData import AIPSUVData as wizAIPSUVData
from Wizardry.AIPSData import AIPSTableRow
from AIPSTask import AIPSTask, AIPSList
from AIPSData import  AIPSImage
import sys, operator, pickle, os

AIPS.userno = 100
fileend = 'PBCOR_IM.fits'

for i in os.listdir('./'):
    if i.endswith(fileend):
        fitld = AIPSTask('FITLD')
        fitld.datain = 'PWD:%s' % i
        fitld.outname = i[:8]
        fitld.go()
