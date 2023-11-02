#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <time.h>
#include <vector>

void load_data(std::string file_name, std::vector<double> *data, std::string delimiter, unsigned int *size = nullptr);
