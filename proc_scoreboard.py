# pre-process all atom pdb
# extract CA in Chain A and save to file
import sys
def main():
	if len(sys.argv) < 2:
		print "Usage python proc_scoreboard.py cluster_file"
		return
	arrSingleVar = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	arrDoubleVar = ['AA','CC','DD','EE','FF','GG','HH','II','KK','LL','MM','NN','PP','QQ','RR','SS','TT','VV','WW','YY']
	arrTripleVar = ['AAA','CCC','DDD','EEE','FFF','GGG','HHH','III','KKK','LLL','MMM','NNN','PPP','QQQ','RRR','SSS','TTT','VVV','WWW','YYY']

	fin = open(sys.argv[1], 'r')
	out_name = sys.argv[1]+'.score'
	fout = open(out_name, 'w')
	for line in fin.readlines():
		line = line.strip()
		strArr = line.split(',')
		groupStr = ''.join(sorted(strArr[2]))
		board = ''
		for i in xrange(0,len(arrTripleVar)):
			if arrTripleVar[i] in groupStr:
				outStr = '0 0 1'
			elif arrTripleVar[i][0:2] in groupStr:
				outStr = '0 1 0'
			elif arrTripleVar[i][0] in groupStr:
				outStr = '1 0 0'
			else:
				outStr = '0 0 0'
			board = ('%s %s') % (board, outStr)

		board = board.lstrip(' ')
		fout.write(('%s,%s,%s,%s\n') % (strArr[0], strArr[2], board, strArr[3]))

	fin.close()
	fout.close()
if __name__=="__main__":
	main()

