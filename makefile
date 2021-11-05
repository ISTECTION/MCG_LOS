SRC	    := src
INCLUDE	:= include
OUTPUT	:= output
CXX 	:= g++
FLAG    := -O2
BOOST   := "D:/Program Files/boost/include/boost-1_77"

RM := del /q /f

compile_running_main: clean compile running

running:
	$(OUTPUT)/main.exe

compile:
	$(CXX) $(FLAG) $(SRC)/main.cpp -o $(OUTPUT)/main.exe -I$(INCLUDE)


compile_generator:
	$(CXX) $(FLAG) $(SRC)/main_generator.cpp -o $(OUTPUT)/main_generator.exe -I$(INCLUDE)

compile_automation:
	$(CXX) $(FLAG) $(SRC)/main_automation.cpp -o $(OUTPUT)/main_automation.exe -I$(INCLUDE)

run_test: clean_test compile_test
	test/output/MCG_none.exe
	test/output/MCG_diag.exe
	test/output/LOS_none.exe
	test/output/LOS_diag.exe

compile_test:
	$(CXX) test/src/MCG_none.cpp -o test/output/MCG_none.exe -I$(INCLUDE) -I$(BOOST)
	$(CXX) test/src/MCG_diag.cpp -o test/output/MCG_diag.exe -I$(INCLUDE) -I$(BOOST)
	$(CXX) test/src/LOS_none.cpp -o test/output/LOS_none.exe -I$(INCLUDE) -I$(BOOST)
	$(CXX) test/src/LOS_diag.cpp -o test/output/LOS_diag.exe -I$(INCLUDE) -I$(BOOST)

clean:
	$(RM) $(OUTPUT)\*.exe

clean_test:
	$(RM) test\output\*.exe
