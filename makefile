FokkerPlanckSolver	:	main.o solver.o IDLinterface.o
				g++ -shared -o FokkerPlanckSolver.so main.o solver.o IDLinterface.o
main.o			:	main.cpp main.h MyTypes.h solver.h
				g++ -c -O3 -fPIC -D LINUX main.cpp
solver.o		:	solver.cpp MyTypes.h solver.h const.h
				g++ -c -O3 -fPIC -D LINUX solver.cpp
IDLinterface.o		:	IDLinterface.cpp main.h
				g++ -c -O3 -fPIC -D LINUX IDLinterface.cpp