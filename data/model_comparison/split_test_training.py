import random

input_file = "Data_preparation/tb_data_10.txt"

train_file = "tb_data_10_train.txt"
test_file = "tb_data_10_test.txt"

random.seed(31)

# import data
with open(input_file, "r") as f:
    lines = [line for line in f if line.strip]

# random permutation of lines
random.shuffle(lines)

# 80-20 split
split = int(0.8* len(lines))

train_lines = lines[:split]
test_lines = lines[split:]

# save
with open(train_file, "w") as f:
    f.writelines(train_lines)

with open(test_file, "w") as f:
    f.writelines(test_lines)