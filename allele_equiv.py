#allele equivalences, check out all the antigen2 
import pandas as pd


equiv_A = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.equiv.A.csv")
equiv_B = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.equiv.B.csv")
equiv_C = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.equiv.C.csv")
equiv_Q = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.equiv.Q.csv")
equiv_R = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.equiv.R.csv")

equiv_dict = {"A":equiv_A,
			"B":equiv_B,
			"C":equiv_C,
			"Q":equiv_Q,
			"R":equiv_R
			}



#equiv_dict.to_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table_equiv.csv") 
def split_char_num(char_num):#split the letter and number
	# input: "QR1"
	# output: ["QR",1]
	idx = 0
	while not char_num[idx].isdigit():
		idx += 1
	# print((char_num[0:idx], char_num[idx:]))
	return [char_num[0:idx], int(char_num[idx:])]


def get_combination(arr, length):#arr=["A","B","C"]是个字符串的dict，找出所有长度为length的所有子串. arr is a character, find out all subset of length = length
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

def get_frequency(char_num, nmdp_5l, ethnics):#[["A",1], ["B",1]] 只考虑一个.sum up the frequency
	selected_data = nmdp_5l
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

def find_equiv(table_equiv,antigen):
	antigen = antigen.replace("DQ","Q")
	antigen = antigen.replace("DR","R")
	antigen = antigen.split(",")
	antigen = [split_char_num(char_num) for char_num in antigen]

	antigen_map = {}
	for char,num in antigen:#把字母下面的数字找到，放到antigen_map字典里面.find out all number under character and put them into antigen_map
		if char not in antigen_map:
			antigen_map[char] = []
		antigen_map[char].append(num)
		print("######",num)
	antigen_char = sorted(antigen_map.keys())  #antigen_mao.keys(): put all dict keys into a list. sorted(), sort the list.
	antigen_char_comb = []
	for i in range(1,len(antigen_char)+1):
		antigen_char_comb += get_combination(antigen_char,i)


	print(antigen)
	print()
	print(antigen_char)
	print()
	print(antigen_char_comb)
	print()
	print(antigen_map)
	print()
	char_num_equiv = []
	char_num_comb = []
	all_antigen = []
	for char_comb in antigen_char_comb:
		char_num_comb += get_char_num_comb(char_comb,antigen_map)
	print("111",char_comb)
	for char_num in char_num_comb:
		print("###",char_num)
	for char,num in antigen:#用list可能更好
		if char == "A":
			print("hahaha,",num)
			print("hahaha,",char)
			dict_equiv_a = equiv_dict["A"]
			equiv_a_bool = dict_equiv_a["Ag"].isin([num])
			equiv_a = dict_equiv_a[equiv_a_bool].dropna(axis=1)#equiv_a is a dataframe
			char_num_equiv.append(equiv_a)
			equiv_a = equiv_a.values.tolist()[0]
			print(equiv_a)
			if len(equiv_a) >= 1:
				equiv_a = ["A" + str(int(x)) for x in equiv_a[2:]]
				all_antigen.extend(equiv_a) 
			else:
				equiv_a = [char + str(num)]
				all_antigen.extend(equiv_a) 

			#equiv_a = ["A" + str(x) for x in equiv_a]
			print(equiv_a)
			#print(all_antigen)
			
			

			#print("equiv_str_a is,", equiv_str_a)

		if char == "B":
			dict_equiv_b = equiv_dict["B"]
			equiv_b_bool = dict_equiv_b["Ag"].isin([num])
			equiv_b = dict_equiv_b[equiv_b_bool].dropna(axis=1)
			char_num_equiv.append(equiv_b)
			equiv_b = equiv_b.values.tolist()[0]
			if len(equiv_b) >= 1:
				equiv_b = ["B" + str(int(x)) for x in equiv_b[2:]]
				all_antigen.extend(equiv_b) 

			else:
				equiv_b = [char + str(num)]
				all_antigen.extend(equiv_b) 
			equiv_str_b = []
			print(equiv_b)

		if char == "C":
			dict_equiv_c = equiv_dict["C"]
			equiv_c_bool = dict_equiv_c["Ag"].isin([num])
			equiv_c = dict_equiv_c[equiv_c_bool].dropna(axis=1)
			char_num_equiv.append(equiv_c)
			equiv_c = equiv_c.values.tolist()[0]

			if len(equiv_c) >= 1:
				equiv_c = ["C" + str(int(x)) for x in equiv_c[2:]]
				all_antigen.extend(equiv_c) 
			else:
				equiv_c = [char + str(num)]
				all_antigen.extend(equiv_c)  
			equiv_str_c = []
			print(equiv_c)

		if char == "Q":
			dict_equiv_q = equiv_dict["Q"]
			equiv_q_bool = dict_equiv_q["Ag"].isin([num])
			equiv_q = dict_equiv_q[equiv_q_bool].dropna(axis=1)
			char_num_equiv.append(equiv_q)
			equiv_q = equiv_q.values.tolist()
			#print(equiv_q)
			if len(equiv_q) >= 1:
				equiv_q = ["Q" + str(int(x)) for x in equiv_q[2:]]
				all_antigen.extend(equiv_q)
			else:
				equiv_q = [char + str(num)]
				all_antigen.extend(equiv_q) 
			
			equiv_str_q = []
			print(equiv_q)

		if char == "R":
			dict_equiv_r = equiv_dict["R"]
			equiv_r_bool = dict_equiv_r["Ag"].isin([num])
			equiv_r = dict_equiv_r[equiv_r_bool].dropna(axis=1)
			char_num_equiv.append(equiv_r)
			equiv_r = equiv_r.values.tolist()[0]
			if len(equiv_r) >= 1:
				equiv_r = ["R" + str(int(x)) for x in equiv_r[2:]]
				all_antigen.extend(equiv_r) 
			else:
				equiv_r = [char + str(num)]
				all_antigen.extend(equiv_r)  
			equiv_str_r = []
			print(equiv_r)
	print()
	print("final_antigen,",all_antigen)
	#sum_list = equiv_a + equiv_b + equiv_c + equiv_q + equiv_r
	#print(sum_list)
	
	#for item in sum_list[0]:


	#	if len(item) > 2 :
	#		item = item[2:]
	#		print("1",item)	
	#	else:
	#		print("2")
	antigen2 = all_antigen
	antigen2 = ','.join(antigen2)
	antigen2 = str(antigen2)
	print(antigen2)
	return antigen2
	

	
def func(nmdp_5l, antigen, ethnics):
	print(antigen)
	antigen = antigen.replace("DQ","Q")
	antigen = antigen.replace("DR","R")
	antigen = antigen.split(",")
	antigen = [split_char_num(char_num) for char_num in antigen]
	print("####")
	print(type(antigen))
	# antigen = [["A",1], ["A",2], ... ,["DR",18]]
	antigen_map = {}
	for char,num in antigen:#把字母下面的数字找到，放到antigen_map字典里面.find out all number under character and put them into antigen_map
		if char not in antigen_map:
			antigen_map[char] = []
		antigen_map[char].append(num)
	antigen_char = sorted(antigen_map.keys())  #antigen_mao.keys(): put all dict keys into a list. sorted(), sort the list.
	antigen_char_comb = []
	for i in range(1,len(antigen_char)+1):
		antigen_char_comb += get_combination(antigen_char,i)
	print(antigen)
	print()
	print(antigen_char)
	print()
	print(antigen_char_comb)
	print()
	print(antigen_map)
	print()
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
		print(get_frequency(char_num, nmdp_5l, ethnics))
		if len(char_num) == 1:
			#print(len(char_num))
			S1 = get_frequency(char_num, nmdp_5l, ethnics).sum(axis=1)
			L1.append(S1)
			print("S1=",S1)
		if len(char_num) ==2:
			#print(len(char_num))
			S2 = get_frequency(char_num, nmdp_5l , ethnics).sum(axis=1)
			L2.append(S2)
			print("S2=",S2)
			
		if len(char_num) ==3:
			#print(len(char_num))
			
			S3 = get_frequency(char_num,nmdp_5l,ethnics).sum(axis=1)
			L3.append(S3)
			print("S3=",S3)
		if len(char_num) ==4:
			#print(len(char_num))
			
			S4 = get_frequency(char_num,nmdp_5l,ethnics).sum(axis=1)
			L4.append(S4)
			print("S4=",S4)
		if len(char_num) ==5:
			#print(len(char_num))
			
			S5 = get_frequency(char_num,nmdp_5l,ethnics).sum(axis=1)
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
	
	antigen_5l = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/table.5l.nmdp.ag.csv")
	ethnics = pd.read_csv("/Users/muwuxu/Documents/Loran_Gragert/cpra_python/data/ethnic.weights_python.csv")
	func(antigen_5l, find_equiv(equiv_dict,"B4,A1,A2,C5,DQ15,DR4"), ethnics)



