import libra_py.packages.cp2k.methods as CP2K_methods

params = {"logfile_name": "step_1.log", "tolerance": 0.00, "number_of_states": 15}

min_ = 1000
max_ = 0
for i in range(1000, 4000):
   print('-------------------------------', i, '-----------------------------')
   params.update({"logfile_name":F"all_logfiles/step_{i}.log"});
   tmp = CP2K_methods.read_cp2k_tddfpt_log_file(params)
   print(tmp)

   
   if (int(tmp[0][0])) < 0 or (tmp[1][0][0][0]) < 2 or (tmp[1][0][0][1]) > 23:
      print('Flago not found', i)
   for j in range(len(tmp[1])):
       for k in range(len(tmp[1][j])):
           min_ = min(tmp[1][j][k][0], min_)
           max_ = max(tmp[1][j][k][1], max_)

   if len(tmp[0]) != 15:
            print("not one")


print('Min index is:', min_)
print('Max index is:', max_)
 
