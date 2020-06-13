SYSTEM:=x86-64_osx

CXX:=clang++

CPLEX_LN_FLAGS:=-lconcert -lilocplex -lcplex
CPLEX_LN_FLAGS+=-m64 -lm -lpthread 
CPLEX_LN_FLAGS+=-framework CoreFoundation 
CPLEX_LN_FLAGS+=-framework IOKit
CPLEX_LN_FLAGS+=-stdlib=libc++

CPLEX_HOME:=/Applications/CPLEX_Studio128
UTIL_HOME:=/Users/ckhan/workspace/util
ROU_HOME:=/Users/ckhan/workspace/BnC_CPLEX/BnC_CPLEX

SRC_DIR:=./src
BUILD_DIR:=./build
LIB_DIR:=./lib

LIB_FORMAT:=static_pic
SRC_EXT:=cpp
OBJ_EXT:=o

CXX_OPT:=-std=c++14 -DIL_STD
#####################################################################
#
# Set paths
#
CPLEX_DIR:=$(CPLEX_HOME)/cplex
CPLEX_INC_DIR:=$(CPLEX_DIR)/include
CPLEX_LIB_DIR:=$(CPLEX_DIR)/lib/$(SYSTEM)/$(LIB_FORMAT)

CONCERT_DIR:=$(CPLEX_HOME)/concert
CONCERT_INC_DIR:=$(CONCERT_DIR)/include
CONCERT_LIB_DIR:=$(CONCERT_DIR)/lib/$(SYSTEM)/$(LIB_FORMAT)

UTIL_NAME:=util
UTIL_INC_DIR:=$(UTIL_HOME)/include
UTIL_LIB_DIR:=$(UTIL_HOME)/lib

ROUTER_NAME:=router
ROU_INC_DIR:=$(ROU_HOME)/include
ROU_LIB_DIR:=$(ROU_HOME)/lib
SEQ_NAME:=sequencer

#####################################################################
MAIN_CPP:=main.$(SRC_EXT)
MAIN_OBJ:=$(BUILD_DIR)/main.$(OBJ_EXT)

#####################################################################
SOURCES:=$(shell find $(SRC_DIR) -type f -name "*.$(SRC_EXT)")
OBJECTS:=$(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.$(SRC_EXT)=.$(OBJ_EXT)))

#####################################################################
CXX_INC_DIRS:=-I$(CPLEX_INC_DIR) -I$(CONCERT_INC_DIR) -I$(UTIL_INC_DIR) -I$(ROU_INC_DIR)
CXX_LIB_DIRS:=-L$(CPLEX_LIB_DIR) -L$(CONCERT_LIB_DIR) -L$(UTIL_LIB_DIR) -L$(ROU_LIB_DIR) 

CXX_FLAGS:=$(CXX_OPT) $(CXX_INC_DIRS)
CXX_LN_FLAGS:=$(CPLEX_LN_FLAGS) -l$(UTIL_NAME) -l$(ROUTER_NAME) -l$(SEQ_NAME)

#####################################################################
PROGRAM:=recommender

#####################################################################

all: $(PROGRAM)

echoTest:
	@echo "echo TEST"
	@echo $(SOURCES)

$(PROGRAM): $(MAIN_OBJ) $(OBJECTS)
	$(CXX) $(CXX_FLAGS) $(CXX_LIB_DIRS) $(MAIN_OBJ) $(OBJECTS) -o $(PROGRAM) $(CXX_LN_FLAGS) 

$(MAIN_OBJ): $(MAIN_CPP)
	@mkdir -p $(dir $@)
	@echo "Start compile"
	$(CXX) -c $(CXX_FLAGS) $< -o $@
	@echo "Build Success!!!"
	@echo ""

#Compile each .cpp file
$(BUILD_DIR)/%.$(OBJ_EXT): $(SRC_DIR)/%.$(SRC_EXT)
	@mkdir -p $(dir $@)
	@echo "Start compile"
	@echo "Target file:" $@ "; Pre-req:" $<
	$(CXX) -c $(CXX_FLAGS) $< -o $@
	@echo "Build Success!!!"
	@echo ""

clean:
	@rm -rf $(BUILD_DIR)
	@rm $(PROGRAM)
