
CXX 		=g++
LIBS		=`root-config --libs`  

CXXFLAGS 	=-g3 -std=c++11  `root-config --cflags` 
INCPATH     =include/

TARGET		= mvaDataMaker
all:	$(TARGET) 

mvaTrainer2 : photonID_bdt_MLP_classification.C 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaTester : testMVAModels.cc MVATrainer.h
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaDataMaker : src/mvaDataMaker.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

mvaTrainer: makeMVAModels.cc  MVATrainer.h MVATrainer.cc 
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@.exe $<  $(LIBS)

clean :
	@rm *.exe
