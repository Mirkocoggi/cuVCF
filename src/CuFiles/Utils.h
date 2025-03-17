#ifndef UTILS_H
#define UTILS_H

#include <sys/wait.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>

const int FLAG = 0;
const int INT = 1;
const int FLOAT = 2;
const int STRING = 3;
const int INT_ALT = 4;
const int FLOAT_ALT = 5;
const int STRING_ALT = 6;
const int STRING_FORMAT = 8;
const int INT_FORMAT = 9;
const int FLOAT_FORMAT = 10;
const int STRING_FORMAT_ALT = 11;
const int INT_FORMAT_ALT = 12;
const int FLOAT_FORMAT_ALT = 13;

//Securely unzip the file
void unzip_gz_file(char* vcf_filename) {
    // Check the extension ".gz"
    if (strcmp(vcf_filename + strlen(vcf_filename) - 3, ".gz") == 0) {
        pid_t pid = fork();
        if (pid == 0) {
            // Child proces
            execlp("gzip", "gzip", "-df", vcf_filename, nullptr);
            // If execlp fails
            cout<< "ERROR: Failed to execute gzip command" << std::endl;
            exit(EXIT_FAILURE);
        } else if (pid > 0) {
            // Parent proces waits for the child proces
            int status;
            waitpid(pid, &status, 0);
            if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
                // Remove ".gz" from the filename
                char* mutable_vcf_filename = const_cast<char*>(vcf_filename);
                mutable_vcf_filename[strlen(vcf_filename) - 3] = '\0';
            } else {
                cout<< "ERROR: cannot unzip file" << std::endl;
            }
        } else {
            // Error in the fork() call
            cout<< "ERROR: Failed to fork process" << std::endl;
        }
    }
}

string get_filename(string path_filename, string &path_to_filename){
    vector<string> line_el;
    path_to_filename = path_filename;
    boost::split(line_el, path_filename, boost::is_any_of("/"));
    return line_el[line_el.size()-1];
}

long get_file_size(string filename){
    return filesystem::file_size(filename);
}




#endif