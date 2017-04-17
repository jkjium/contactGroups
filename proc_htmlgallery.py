import sys

# convert cathpair aligned figures into html table
# cath aligned pair file name: p.2dpmA02.1q0sA02.seq.cb2.png
# given list format: 2dpmA02 1q0sA02

def main():
	if len(sys.argv) < 2:
		print 'Usage: python proc_htmlgallery.py pairlist.txt'
		return

	prefix = ['<html>' ,
					'<head>' ,
					'<link type="text/css" rel="stylesheet" href="style.css">',
					'</head>',
					'<body>',
						'<div id="content">',
							'<table>\n'
			]

	suffix = [
					'\n</table>',
				'</div>',
			'</body>',
		'</html>'
	]

	# table content
	'''
				<tr>
					<td><img src="png/p.1kmoA01.3efmA01.seq.b62.png"></td>
					<td><img src="png/p.1kmoA01.3efmA01.seq.cb2.png"></td>
				</tr>
				<tr>
					<td>1kmoA01 , 3efmA01 - BLOSUM62</td>
					<td>1kmoA01 , 3efmA01 - SCSC</td>
				</tr>
	'''

	pnglistfile = sys.argv[1]
	outfile = pnglistfile+'.html'

	content = []
	with open(pnglistfile) as fp:
		for line in fp:
			line = line.strip()
			titlestr = line.split(' ')
			pngtr = '<tr><td><img src="p.%s.%s.seq.b62.png"></td>' % (titlestr[0], titlestr[1]) + \
					'<td><img src="p.%s.%s.seq.cb2.png"></td></tr>' % (titlestr[0], titlestr[1])
			captr = '<tr><td>%s , %s - BLOSUM62</td>' % (titlestr[0], titlestr[1]) + \
					'<td>%s, %s - SCSC</td></tr>' % (titlestr[0], titlestr[1])
			content.append(pngtr)
			content.append(captr)

	fout = open(outfile, 'w')
	fout.write('\n'.join(prefix))
	fout.write('\n'.join(content))
	fout.write('\n'.join(suffix))
	fout.close()

if __name__ == '__main__':
	main()