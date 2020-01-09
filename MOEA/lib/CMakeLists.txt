# project
#-------------------------------
project(MOEAlib)

# Set public header list (add your headers and source files here))
#-------------------------------
set(${PROJECT_NAME}_headers
                        Chromosome.h
                        Crossover.h
                        Data.h
                        eoETTPEval.h
                        eoEvolutionOperator.h
                        eoSFLA.h
                        eoSFLAInitializer.h
                        ETTPInit.h
                        ETTPKempeChainHeuristic.h
                        ETTPneighbor.h
                        ETTPneighborEval.h
                        ETTPneighborhood.h
                        Exam.h
                        GraphColouringHeuristics.h
                        ITC07TestSet.h
                        Matrix.h
                        MOEA.h
                        moeoETTPEval.h
                        moeoLocalSearchUpdater.h
                        moeoSFLA.h
                        moLASAexplorer.h
                        moSHC.h
                        moSHCexplorer.h
                        Mutation.h
                        Period.h
                        Repairing.h
                        TestSet.h
                        TestSetDescription.h
                        TorontoTestSet.h
                        VertexPriorityQueue.h
)
set(${PROJECT_NAME}_sources
                        Chromosome.cpp
                        Crossover.cpp
                        ETTPKempeChainHeuristic.cpp
                        GraphColouringHeuristics.cpp
                        ITC07TestSet.cpp
                        MainApp.cpp
                        Matrix.cpp
                        MOEA.cpp
                        Mutation.cpp
                        TestSetDescription.cpp
                        TorontoTestSet.cpp
                        VertexPriorityQueue.cpp
)


add_library(MOEAlib SHARED ${${PROJECT_NAME}_headers} ${${PROJECT_NAME}_sources})


# Include ParadisEO libs and Boost Regex lib
target_link_libraries(${PROJECT_NAME} boost_regex-gcc-1_52 eo es moeo cma eoutils ga)

