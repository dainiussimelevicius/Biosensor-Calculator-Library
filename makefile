calculator_tester: calculator_tester.o libbiosensor_calculator.so
	gcc calculator_tester.o -Wl,-rpath,. -L. -lbiosensor_calculator -o calculator_tester

libbiosensor_calculator.so: biosensor_calculator.o explicit_calculator.o implicit_calculator.o utils.o
	gcc -shared biosensor_calculator.o explicit_calculator.o implicit_calculator.o utils.o -o libbiosensor_calculator.so

calculator_tester.o: calculator_tester.c biosensor_calculator.h biosensor_information.h
	gcc -c calculator_tester.c -o calculator_tester.o

biosensor_calculator.o: biosensor_calculator.c biosensor_calculator.h biosensor_information.h
	gcc -c -fPIC biosensor_calculator.c -o biosensor_calculator.o

explicit_calculator.o: explicit_calculator.c explicit_calculator.h biosensor_information.h utils.h constants.h
	gcc -c -fPIC explicit_calculator.c -o explicit_calculator.o

implicit_calculator.o: implicit_calculator.c implicit_calculator.h biosensor_information.h utils.h constants.h
	gcc -c -fPIC implicit_calculator.c -o implicit_calculator.o

utils.o: utils.c utils.h
	gcc -c -fPIC utils.c -o utils.o

clean:
	rm libbiosensor_calculator.so biosensor_calculator.o explicit_calculator.o implicit_calculator.o utils.o calculator_tester.o calculator_tester
