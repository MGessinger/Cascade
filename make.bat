gcc -Wall -pedantic -Wextra -o cascade.exe cascade.c acb_ode.c juliaInterface.c -lflint -larb -lgmp
cascade.exe data/testing.txt 10 50 0.1
pause
