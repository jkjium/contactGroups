import sys
import json

def main():
	if len(sys.argv) < 2:
		print 'Usage: python proc_pfamscan_json.py 1bjp_A.json'
		return

	with open(sys.argv[1]) as fp:
		line = fp.readline()

	jsonlist = json.loads(line.strip())

	for i in xrange(0, len(jsonlist)):
		j = jsonlist[i]

		name = j['seq']['name']
		start = j['seq']['from']
		end = j['seq']['to']

		pfam = j['acc'][0:7]

		seqstr = j['align'][3]
		align = seqstr.split()[1]

		seq = align.replace('-', '').replace('.', '')

		print '%s %s %s %s %s' % (name, start, end, pfam, seq)

if __name__ == '__main__':
	main()