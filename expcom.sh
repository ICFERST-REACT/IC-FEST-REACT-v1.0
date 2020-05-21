g++ ../phreeqc/Mjd5.cpp -I ../include/ -L ../lib/ -l phreeqcrm
export LD_LIBRARY_PATH=../lib
./a.out
