#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

vector<vector<double>> GetTable(string filename, vector<string> columns) {
  string line;
  ifstream file(filename);

  vector<int> labels;
  map<int, int> colmap;

  int counter = 0;
  int numberoflines = 1;

  // check if file is open
  if (file.is_open()) {
    // printf("File is open\n");

    // read number of lines
    while (!file.eof()) {
      getline(file, line);
      ++numberoflines;
    }
    // cout << numberoflines << endl;
    // printf("Number of lines %i\n", numberoflines);
    file.clear();
    file.seekg(0, ios::beg);
    // getline(file, line);
    // cout << line << endl;
    vector<vector<double>> output(numberoflines - 50,
                                  vector<double>(columns.size()));

    // read lines until eof
    while (!file.eof()) {
      // increase line counter
      ++counter;
      // cout << counter << endl;

      // read a line
      getline(file, line);

      // check if eof
      if (file.eof())
        break;

      // if line is 47 read the column names
      if (counter == 47) {
        // load the current line as stream
        istringstream iss(line);

        // split the line and save in vector
        // vector<string> labels;
        int colcounter = 0;
        // cout << "Col idx: ";
        do {
          string sub;
          iss >> sub;
          vector<string>::iterator it =
              find(columns.begin(), columns.end(), sub);
          if (it != columns.end()) {
            colmap[colcounter - 1] = it - columns.begin();
            labels.push_back(colcounter - 1);
            // cout << colcounter << " ";
          }
          ++colcounter;
        } while (iss);
        // cout << endl;
      }

      if (counter > 48) {
        istringstream iss(line);
        int colcounter = 0;
        int seccolcounter = 0;
        do {
          string sub;
          iss >> sub;
          vector<int>::iterator iti =
              find(labels.begin(), labels.end(), colcounter);
          if (iti != labels.end()) {
            // cout << counter << " " << counter - 49 << " " << sub;
            // output[counter - 49].push_back(stod(sub));
            output[counter - 49][seccolcounter] = stod(sub);
            ++seccolcounter;
          }
          ++colcounter;
        } while (iss);
        // cout << endl;
      }
    }
    return output;
  }
  vector<vector<double>> output(1, vector<double>(1));
  output[0][0] = 0.0;
  return output;
}
