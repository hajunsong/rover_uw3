#include "fileio.h"

void load_data(std::string file_name, std::vector<double> *data, std::string delimiter, unsigned int *size)
{
	data->clear();
	FILE *fp_in;
	const int buffer = 1000000;
	char *ptr, basic[buffer];
	fp_in = fopen(file_name.c_str(), "r");
	unsigned int row = 0, col = 0;
	while (fgets(basic, buffer, fp_in) != NULL)
	{
		row++;
		ptr = strtok(basic, delimiter.c_str());
		while (ptr != NULL) {
			data->push_back(atof(ptr));
			ptr = strtok(NULL, delimiter.c_str());
			col++;
		}
	}
	fclose(fp_in);
	if(size != NULL){
		size[0] = row;
		size[1] = (col/row);
	}
}
