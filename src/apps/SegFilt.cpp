/*
 * SegFilt.cpp
 * 
 * Copyright 2015 Robert Bakaric <rbakaric@irb.hr>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */




#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<boost/program_options.hpp>

#include<SEG.hpp>


namespace po = boost::program_options;
using namespace std;



void PrintAcknowledgements(){

cout <<"***************************************************************************************"<<"\n";
cout <<"                                 SegFilt - SEG Filter                                  "<<"\n";
cout <<"                                          by                                           "<<"\n";
cout <<"                                    Robert Bakaric                                     "<<"\n\n";
cout <<"CONTACT:                                                                               "<<"\n";
cout <<" Code written and maintained by Robert Bakaric,                                        "<<"\n";
cout <<" email: rbakaric@irb.hr , bakaric@evolbio.mpg.de                                       "<<"\n\n";
cout <<"ACKNOWLEDGEMENT:                                                                       "<<"\n";
cout <<"   Wootton, J. C. and S. Federhen (1993).  Statistics of local complexity in amino     "<<"\n";
cout <<"   acid sequences and sequence databases.  Computers and Chemistry 17:149-163.         "<<"\n";
cout <<"                                                                                       "<<"\n";
cout <<"LICENSE:                                                                               "<<"\n";
cout <<" The program is distributed under the GNU General Public License. You should have      "<<"\n";
cout <<" received a copy of the licence together  with this software. If not, see              "<<"\n";
cout <<" http://www.gnu.org/licenses/                                                          "<<"\n";
cout <<"***************************************************************************************"<<"\n\n\n";
}


template <typename INT, typename CHARA>
po::variables_map SetOptions(INT& argc, CHARA& argv){

    try {
        int opt;
        string version = "1.0";
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version,v", "print version information")
            ("input-file,i", po::value< string >(), "Fasta input file")
            ("window,w", po::value< int >()->default_value(12), "SEG window size")
            ("locut,L", po::value< double >()->default_value(2.2), "Low complexity cutoff (starter)")
            ("hicut,H", po::value< double >()->default_value(2.5), "High complexity cutoff (starter)")
            ("maxXes,x", po::value< int >()->default_value(0), "Maximum number of xxx  symboles tolerated (dynamically defined if left unchanged)")
            ("maxTrim,t", po::value< int >()->default_value(100),"Maximum trimming of raw segment")
        ;

        po::positional_options_description p;
        p.add("input-file,i", -1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);
    
        if (vm.count("help")) {
            cout << "Usage: ./program [options]\n\n";
            PrintAcknowledgements();
            cout << desc;
            exit(0);
        }else if (vm.count("version")) {
            cout << "Program version:  " << version << "\n";
            exit(0);
        }
        if (!vm.count("input-file")){
            cout << "Input file is not defined \n";
            exit(0);
        }
        return vm;    
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        exit(0);
    }    
}

template <typename Tseg>
void FilterFile(const string file, Tseg& seg){
  
  fstream fs;
  
  fs.open (file.c_str(), ios::in);
  
  if ( !fs.is_open())
    throw  runtime_error ("Cannot open file: " + file );

  string str;
  char c;
  int i=0, f=0;
  
  while( fs.good()) {
    c = fs.get();
    if(fs.eof())
       break;
    if (c == '>'){
		if (str.size()> 0){
            str.resize(i);
            string s = seg.SegFilt(str);
            cout << s << endl;
            str.resize(0);
            i=0;
		}
        f=1;
        cout<<c;
     }else if( f == 1 && c == '\n'){
        f=0;
        cout<<c;
     }else if( f == 1 ){
        cout<<c;
     }else if (f == 0 && c !='\n'){
		 if(i == str.size() )
		     str.resize((str.size()+1)*2);
		 str[i++] = c;
	 }
  }

/* Flush */
str.resize(i);
string s = seg.SegFilt(str);
cout << s << endl;
fs.close();
}



int main (int argc, char** argv){

/* Get cmd Options */  
  po::variables_map arg;
  arg = SetOptions(argc,argv);
  

/* Parse Options */
  unordered_map<string, int> args;
  
  args["window"] = static_cast<int>(arg["window"].as<int>());
  args["hicut"] =  static_cast<double>(arg["hicut"].as<double>());
  args["locut"] =  static_cast<double>(arg["locut"].as<double>());
  args["maxtrim"] =  static_cast<int>(arg["maxTrim"].as<int>());
  args["maxXes"] =  static_cast<int>(arg["maxXes"].as<int>());

  try{

/* Make object */
    SEG<int> seg(args);
  
/* Filter Sequences */  
    FilterFile(arg["input-file"].as<string>(), seg);
    
  }catch(runtime_error& e){
    cerr << e.what() << "\n";
  }
  
  return 0;
}
