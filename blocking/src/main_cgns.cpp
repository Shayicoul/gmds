//
// Created by calderans on 30/06/23.
//
#include <gmds/blocking/CGNSWriter3D.h>

using namespace gmds;
using namespace blocking;
/*----------------------------------------------------------------------------*/
int main(int argc, char** argv)
{
	std::cout << "CGNS writer:" << std::endl;
	if (argc != 4 )
	{
		std::cout << "Usage: CGNSWriter <input_magix_mesh>.vtk <path_for_output> <output_name>.cgns" << std::endl;
		exit(0);
	}
	std::string param_file(argv[1]);
	std::string path(argv[2]);
	std::string filename(argv[3]);

	CGNSWriter3D writer3D;
	writer3D.write(param_file,filename,path);

}