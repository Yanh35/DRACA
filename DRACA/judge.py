import sys
import linecache
dic_breast={}
dic_lung={}
dic_colorectal={}
for i in range(1,1680,1):
    lncrna=linecache.getline('../data/lncrna-id.txt',i)
    score_breast=linecache.getline(sys.argv[2],i)
    dic_breast[lncrna.strip()]=score_breast.strip()
    score_lung=linecache.getline(sys.argv[3],i)
    dic_lung[lncrna.strip()]=score_lung.strip()
    score_colorectal=linecache.getline(sys.argv[4],i)
    dic_colorectal[lncrna.strip()]=score_colorectal.strip()
if float(dic_breast[sys.argv[1]]) >= 0.785:
    print sys.argv[1],'is related to breast cancer'
else:
    print sys.argv[1],'is not related to breast cancer'
if float(dic_lung[sys.argv[1]]) >= 0.965:
    print sys.argv[1],'is related to lung cancer'
else:
    print sys.argv[1],'is not related to lung cancer'
if float(dic_colorectal[sys.argv[1]]) >= 0.815:
    print sys.argv[1],'is related to colorectal cancer'
else:
    print sys.argv[1],'is not related to colorectal cancer'
