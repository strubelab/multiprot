import re

n_chains = 2
filenames = []

with open("filenames.txt") as f:
	for line in f:
		filenames.append(line.rstrip('\n'))

# Count residues
res_count = 0
with open(filenames[0]) as f:
	for line in f:
		if re.search(r' CA ', line):
			res_count += 1

print 'Residue count:', res_count

nchains = int(raw_input("Number of chains: "))
chainsl = []	# list for chain lengths
for i in range(nchains):
	chainsl.append(int(raw_input('Length of chain ' + str(i+1) + ': ')))

print ''

for filename in filenames:

	data = []

	# Copy data from file
	with open(filename) as f:
		for line in f:
			data.append(line.rstrip('\n'))

	##print filename
	##print 'ORIGINAL'
	##for d in data[:20]:
	##	print d
	##print ''
	
	chains = []		# list for each chain atom entries
	hetatm = []		# list for each chain hetatm entries
	res_chain = []	# list for last residue in each chain
	chain_identif = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']	# list for chain identifiers

	ncount = 0
	index = 0		# index for the entry number

	for i in range(nchains):
		# chainsl.append(int(raw_input('Length of chain ' + str(i+1) + ': ')))
		chains.append([])
		hetatm.append([])
		res_chain.append([0])

		chainl = chainsl[i]		# length of the current chain
		cacount = 0		# counts the residues

		# Storing the chains in separate elements of 'chains' list
		for j in range(index,len(data)):
			index += 1
			
			char = re.search(r'[A-Z][A-Z][A-Z] [AB] ', data[j]).group()
			# if NAP, add to hetatm list with the appropriate chain identifier
			if re.search(r'NAP', data[j]):
				hetatm[i].append(data[j].replace(char, char[0:4]+chain_identif[i]+' '))
			# if OXT and at the end of the chain, add to the chain list
			elif re.search(r' OXT ', data[j]) and (cacount == chainl):
				chains[i].append(data[j].replace(char, char[0:4]+chain_identif[i]+' '))
			# if not OXT, add to the chain list
			elif not re.search(r' OXT ', data[j]):
				chains[i].append(data[j].replace(char, char[0:4]+chain_identif[i]+' '))

			residue = data[j][17:20]
			if re.search(r' CA ', data[j]):
				cacount += 1
			if cacount == chainl:
				# if it is the last entry, break ... not sure about this one
				if j == len(data)-1:
					break
				# if the next residue is different, and it is not a NAP atom
				elif residue != data[j+1][17:20] and data[j+1][17:20] != 'NAP':
					break 

	'''
	for i in range(nchains):
		for d in chains[i][:20]:
			print d
	'''

	n = 1		# initialising entry number

	for i in range(nchains):		
		# Replace entry number and residue sequence number
		res_old = chains[i][0][23:26]	# initialize to the first residue number (should be 1 in most cases)
		res_new = 1
		for j in range(len(chains[i])):

			entry_num = re.search(r'ATOM\s+\d+', chains[i][j])
			entry_res = re.search(r'[A-Z][A-Z][A-Z]\s[AB]\s+' + res_old, chains[i][j])
			
			if entry_res:
				# replace residue number
				chains[i][j] = chains[i][j].replace(entry_res.group(), entry_res.group()[0:5] + str(res_new).rjust(4))
				# replace entry number
				chains[i][j] = chains[i][j].replace(entry_num.group(), entry_num.group()[0:5] + str(n).rjust(6))
				n += 1
				continue
			
			# If the residue number changed (new residue)
			res_old = chains[i][j][23:26]
			res_new += 1
			entry_res = re.search(r'[A-Z][A-Z][A-Z]\s[AB]\s+' + res_old, chains[i][j])
			if entry_res:
				# replace residue number
				chains[i][j] = chains[i][j].replace(entry_res.group(), entry_res.group()[0:5] + str(res_new).rjust(4))
				# replace entry number
				chains[i][j] = chains[i][j].replace(entry_num.group(), entry_num.group()[0:5] + str(n).rjust(6))
				n += 1
			else:
				raise AttributeError(str(n) + " : chain " + str(i) + ", entry " + str(j) + ", res_old: "+res_old)

		# Add TER entry at the end of the chain
		# FIX FOR CHAINS WITH MORE THAN 1000 RESIDUES
		aa = re.search(r'[A-Z][A-Z][A-Z]\s[AB]\s+', chains[i][-1])
		chains[i].append('TER    ' + str(n) + aa.group()[0:3].rjust(9) + aa.group()[4].center(3) + str(res_new))
		n += 1

		# last residue number in the current chain, plus 2 to add the hetatm entries
		res_chain[i] = res_new+2
	
	# Join the chains in a single list
	complete = []
	for i in range(nchains):
		complete += chains[i]

	# Add HETATOM entries
	for i in range(nchains):
		res_new = res_chain[i]
		for j in range(len(hetatm[i])):
			entry_num = re.search(r'ATOM\s+\d+', hetatm[i][j])
			entry_res = re.search(r'NAP [AB]\s+\d+', hetatm[i][j])

			# replace residue number and entry number
			complete.append(hetatm[i][j].replace(entry_res.group(), entry_res.group()[0:5] + str(res_new).rjust(4)).replace(entry_num.group(), 'HETATM' + str(n).rjust(5)))
			n += 1

	##print ''
	##print 'FIXED'
	##for d in complete[:20]:
	##	print d
	##print "..."
	##for d in complete[-20:]:
	##	print d
	##print ''

	write_name = filename[:-4]+"_fixed.pdb"
	with open(write_name,"w") as f:
		for d in complete:
			f.write(d + '\n')
