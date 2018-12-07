#allele equivalences, check out all the antigen 
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
	for char_comb in antigen_char_comb:
		char_num_comb += get_char_num_comb(char_comb,antigen_map)
	for char_num in char_num_comb:
		print(char_num)
	for char,num in char_num:#用list可能更好
		if char == "A":
			dict_equiv_a = equiv_dict["A"]
			equiv_a_bool = dict_equiv_a["Ag"].isin([num])
			equiv_a = dict_equiv_a[equiv_a_bool].dropna(axis=1)
			char_num_equiv.append(equiv_a)
			print(equiv_a)
			
		if char == "B":
			dict_equiv_b = equiv_dict["B"]
			equiv_b_bool = dict_equiv_b["Ag"].isin([num])
			equiv_b = dict_equiv_b[equiv_b_bool].dropna(axis=1)
			char_num_equiv.append(equiv_b)
			print(equiv_b)

		if char == "C":
			dict_equiv_c = equiv_dict["C"]
			equiv_c_bool = dict_equiv_c["Ag"].isin([num])
			equiv_c = dict_equiv_c[equiv_c_bool].dropna(axis=1)
			char_num_equiv.append(equiv_c)
			print(equiv_c)

		if char == "Q":
			dict_equiv_q = equiv_dict["Q"]
			equiv_q_bool = dict_equiv_q["Ag"].isin([num])
			equiv_q = dict_equiv_q[equiv_q_bool].dropna(axis=1)
			char_num_equiv.append(equiv_q)
			print(equiv_q)

		if char == "R":
			dict_equiv_r = equiv_dict["R"]
			equiv_r_bool = dict_equiv_r["Ag"].isin([num])
			equiv_r = dict_equiv_r[equiv_r_bool].dropna(axis=1)
			char_num_equiv.append(equiv_r)
			print(equiv_r)



	print("char_num_equiv:",char_num_equiv)


if __name__=="__main__":
	find_equiv(equiv_dict,"B4,A2,C5,DQ15,DR4")



