gcc -O2 -Wall -pedantic -Wextra -s -o cascade.exe cascade.c acb_ode.c -lflint -larb -lgmp
monodromy.exe data/testing.txt 10 50 0.1
pause
