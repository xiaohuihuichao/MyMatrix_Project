CMAKE_MINIMUM_REQUIRED( VERSION 3.0 )



PROJECT( Linear_Regression )




AUX_SOURCE_DIRECTORY( src DIR_LinearRegression )
AUX_SOURCE_DIRECTORY( . DIR_LinearRegression_MAIN )

INCLUDE_DIRECTORIES( include ../../MyMatrix/include)

FILE(GLOB_RECURSE LinearRegression_HEADERS *.h *.hpp)
SOURCE_GROUP("Include" FILES ${LinearRegression_HEADERS} ) 



ADD_EXECUTABLE( 01.Linear_Regression ${DIR_LinearRegression_MAIN}  
				${DIR_LinearRegression} 
				${LinearRegression_HEADERS} )
