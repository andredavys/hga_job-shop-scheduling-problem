
file_instances = open('instances.txt', 'r')
a = 0
while True:
	line = file_instances.readline()
	line_split = line.split()
	if len(line_split) > 0:
		if line_split[0] == 'instance':
			name_file = line_split[1] + '_instance.jsp'
			new_file  = open(name_file, 'w')
			for i in range(4):
				line = file_instances.readline()
			line_split = line.split()
			count_lines = line_split[0]
			new_file.write(line)
			for i in range(int(count_lines)):
				new_file.write(file_instances.readline())n
				
			print('file', name_file, 'created')
		if line_split[0] == 'TERMINOU':
			break