import pandas as pd


#os.chdir('/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data')

#csv = os.listdir()

# antigen = A1, A2 ,B1 ,C2, Q1,R1

#path = '/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/'
#	filename = path + file
#   print(filename) 
#   open(filename)
#input = "A1,A2,A23,A25,A26,A29,A3,A30,A31,A32,A33,A34,A36,A43,A66,A68,A69,A74,A80,B13,B18,B27,B37,B38,B39,B41,B42,B44,B45,B47,B48,B49,B50,B52,B54,B55,B56,B57,B58,B59,B60,B61,B62,B63,B64,B65,B67,B7,B71,B72,B73,B75,B76,B77,B8,B81,B82,C15,C17,DQ6,DR1,DR10,DR103,DR13,DR14,DR16,                                    DR17,DR18"

def split_char_num(char_num):#split the letter and number
	# input: "QR1"
	# output: ["QR",1]
	idx = 0
	while not char_num[idx].isdigit():
		idx += 1
	# print((char_num[0:idx], char_num[idx:]))
	return [char_num[0:idx], int(char_num[idx:])]

def get_combination(arr, length):#arr=["A","B","C"]是个字符串，找出所有长度为length的所有子串. arr is a character, find out all subset of length = length
	if len(arr)<length:
		return None
	if length==1:
		return [[arr[i]] for i in range(len(arr))]
	combination = []
	for i in range(len(arr)):
		cur = [arr[i]]
		suffix = get_combination(arr[i+1:], length-1)
		if suffix is not None:
			combination += [cur+suffix_comb for suffix_comb in suffix]
	return combination

def get_char_num_comb(arr,antigen_map):    #对于各种排列的情况，eg: AB, A:1,2. B:2,3.能有多少种排列组合. find out all combination of A1,2;B2,3.
	char_num_comb = [[]]
	for char in arr:
		new_char_num_comb = []
		num_list = antigen_map[char]
		for num in num_list:
			for char_num in char_num_comb:
				new_char_num_comb.append(char_num+[[char,num]])
		char_num_comb = new_char_num_comb
	return char_num_comb

def get_frequency(char_num, unos_5l, ethnics):#[["A",1], ["B",1]] 只考虑一个.sum up the frequency
	selected_data = unos_5l
	df_ethnics = pd.DataFrame(ethnics)#input ethnics frequencies
	for char,num in char_num:
		selected_data = selected_data.loc[selected_data[char]==num]
		df_freq = pd.DataFrame( {'CAU': [selected_data.loc[:,'CAU'].sum()],
							    'AFA': [selected_data.loc[:,'AFA'].sum()],
							    'HIS': [selected_data.loc[:,'HIS'].sum()],
							    'API': [selected_data.loc[:,'API'].sum()]})
	print(df_ethnics)
	print(df_freq)
	df_weighted_freq=df_freq.mul(df_ethnics)#multiply antigen frequencies with ethnics frequencies
	return df_weighted_freq

def func(unos_5l, antigen, ethnics):
	antigen = antigen.split(",")
	antigen = [split_char_num(char_num) for char_num in antigen]
	# antigen = [["A",1], ["A",2], ... ,["DR",18]]
	antigen_map = {}
	for char,num in antigen:#把字母下面的数字找到，放到antigen_map字典里面.find out all number under character and put them into antigen_map
		if char not in antigen_map:
			antigen_map[char] = []
		antigen_map[char].append(num)
	antigen_char = sorted(antigen_map.keys())
	antigen_char_comb = []
	for i in range(1,len(antigen_char)+1):
		antigen_char_comb += get_combination(antigen_char,i)
	#print(antigen_char_comb)
	#print()
	#print(antigen_map)
	#print()
	L1 = []
	L2 = []
	L3 = []
	L4 = []
	L5 = []
	char_num_comb = []
	for char_comb in antigen_char_comb:
		char_num_comb += get_char_num_comb(char_comb,antigen_map)
	for char_num in char_num_comb:
		print(char_num)
		print(get_frequency(char_num, unos_5l, ethnics))
		if len(char_num) == 1:
			#print(len(char_num))
			S1 = get_frequency(char_num, unos_5l, ethnics).sum(axis=1)
			L1.append(S1)
			print("S1=",S1)
		if len(char_num) ==2:
			#print(len(char_num))
			S2 = get_frequency(char_num, unos_5l , ethnics).sum(axis=1)
			L2.append(S2)
			print("S2=",S2)
			
		if len(char_num) ==3:
			#print(len(char_num))
			
			S3 = get_frequency(char_num,unos_5l,ethnics).sum(axis=1)
			L3.append(S3)
			print("S3=",S3)
		if len(char_num) ==4:
			#print(len(char_num))
			
			S4 = get_frequency(char_num,unos_5l,ethnics).sum(axis=1)
			L4.append(S4)
			print("S4=",S4)
		if len(char_num) ==5:
			#print(len(char_num))
			
			S5 = get_frequency(char_num,unos_5l,ethnics).sum(axis=1)
			L5.append(S5)
			print("S5=",S5)	

	T1 = sum(L1)
	print("T1:",T1)
	T2 = sum(L2)
	print("T2:",T2)
	T3 = sum(L3)
	print("T3:",T3)
	T4 = sum(L4)
	print("T4:",T4)
	T5 = sum(L5)
	print("T5:",T5)

	CPRA = (1-(1-T1+T2-T3+T4-T5)**2)*100
	print("CPRA:",CPRA)



if __name__=="__main__":
	unos_5l = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.5l.unos.csv")
	ethnics = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/ethnic.weights_python.csv")
	func(unos_5l, "A2", ethnics)




