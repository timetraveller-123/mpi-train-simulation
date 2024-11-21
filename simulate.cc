#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>
#include <utility>
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <deque>
#include "platform_load_time_gen.hpp"
#include <algorithm>
#include <iostream>

using namespace std;
using adjacency_matrix = vector<std::vector<size_t>>;




void insertionSort(int arr[][3]) {
    for (int i = 1; i < 9; ++i) {
        int key1 = arr[i][0];  
        int key2 = arr[i][1];  
        int key3 = arr[i][2];
        int j = i - 1;
        
        while (j >= 0 && arr[j][0] > key1) {
            arr[j + 1][0] = arr[j][0];
            arr[j + 1][1] = arr[j][1];
            arr[j + 1][2] = arr[j][2];
            j = j - 1;
        }

        arr[j + 1][0] = key1;
        arr[j + 1][1] = key2;
        arr[j + 1][2] = key3;
    }
}

void simulate(size_t num_stations, const vector<string> &station_names, const std::vector<size_t> &popularities,
              const adjacency_matrix &mat, const unordered_map<char, vector<string>> &station_lines, size_t ticks,
              const unordered_map<char, size_t> num_trains, size_t num_ticks_to_print, size_t mpi_rank,
              size_t total_processes) {
    
    unordered_map<char, size_t> Num_trains = num_trains;

    const int station_max = 17576;

    vector<int> neighbours[station_max], orders[station_max];
    unordered_map<string, int> stationId;
    char lines[3] = {'g', 'y', 'b'};
    int id = 0;

    for(size_t i = 0; i < num_stations; i++){
        stationId[station_names[i]] = i;
    }

    int buffer[station_max][9][3]; // to hold trains before putting in queue                id, line, station_index in line                             negative index for opposite direction of travel
    int links[station_max][6][4]; //to know trains in links                                 id, line, station index in line, time left in link          negative index for opposite direction of travel
    int platforms[station_max][6][4];  // to know trains in platform                        id, line, station_index in line, time left in platform      negative index for opposite direction of travel
    deque<tuple<int,int,int>> queues[station_max][6]; // queues for each platform           id, line, station_index in line                             negative index for opposite direction of travel


    vector<PlatformLoadTimeGen> gens[station_max];


    //calculate neighbours and orders
    for(size_t i = mpi_rank; i < num_stations; i+=total_processes){
        for(size_t j = 0; j < num_stations; j++){
            if(i == j || mat[i][j] <= 0)continue;
            
            neighbours[i].push_back(j);
            int temp = 0;
            for(size_t k = 0; k < i; k++){
                if(mat[j][k] >0)temp++;
            }
            orders[i].push_back(temp);
            
        }
    }

    // calculate all edges
    vector<pair<int,int>> edges;
    for(size_t i = 0; i < num_stations;i++){
        for(size_t j = 0; j < num_stations; j++){
            if(i == j || mat[i][j] <= 0 || i%total_processes == j%total_processes)continue;
            edges.push_back(make_pair(i,j));
        }
    }
    

    //initialize everything
    for(size_t i = mpi_rank; i < num_stations; i+=total_processes){
        for(size_t j = 0; j < 9; j++){
            
            buffer[i][j][0] = buffer[i][j][1] = buffer[i][j][2] = -1;
        }

        for(int j = 0; j < 6; j++){
            platforms[i][j][0] = platforms[i][j][1] = platforms[i][j][2] = platforms[i][j][3] = -1;
            links[i][j][0] = links[i][j][1] = links[i][j][2] = links[i][j][3] = -1;
            gens[i].push_back(PlatformLoadTimeGen(popularities[i])); 
        }
    }

    
    //begin ticks simulation
    for(size_t t = 0; t < ticks; t++){

        // send and receive link updates for stations within the same process
        for(size_t i = mpi_rank; i < num_stations; i+= total_processes){
            for(size_t j = 0; j < neighbours[i].size(); j++){

                int n = neighbours[i][j];
                size_t rank = n%total_processes;
                if(rank != i%total_processes)continue;

                int msg[3] = {-1,-1, -1};
                if(links[i][j][3] > 0){
                    links[i][j][3] -= 1;
                }
                if(links[i][j][3] == 0){
                    msg[0] = links[i][j][0];
                    msg[1] = links[i][j][1];
                    msg[2] = links[i][j][2];

                    links[i][j][0] = links[i][j][1] = links[i][j][2] = links[i][j][3] = -1;
                }

                
                if(rank == i%total_processes){
                    buffer[n][orders[i][j]][0] = msg[0];
                    buffer[n][orders[i][j]][1] = msg[1];
                    buffer[n][orders[i][j]][2] = msg[2];
                }
            }
        }

        // send and receive link update between two processes
        vector<MPI_Request> requests;
        for(int l = 0; l < (int)(edges.size()); l+=100){
            for(int m = l; m < l+100; m++){
                if(m >= (int)(edges.size()))break;
                auto p = edges[m];
                int source = p.first;
                int dest = p.second;
                if(mpi_rank == source%total_processes){
                    int msg[] = {-1,-1,-1};
                    int j=-1;
                    for(int k = 0; k < (int)(neighbours[source].size()); k++){
                        if(neighbours[source][k] == dest){
                            j=k;
                            break;
                        }
                    }

                    if(links[source][j][3] > 0){
                        links[source][j][3] -= 1;
                    }
                    if(links[source][j][3] == 0){
                        msg[0] = links[source][j][0];
                        msg[1] = links[source][j][1];
                        msg[2] = links[source][j][2];

                        links[source][j][0] = links[source][j][1] = links[source][j][2] = links[source][j][3] = -1;
                        
                    }

                    MPI_Request r;
                    MPI_Isend(msg, 3, MPI_INT,dest%total_processes, source*num_stations + dest, MPI_COMM_WORLD, &r); //send
                    requests.push_back(r);


                }else if(mpi_rank == dest%total_processes){
                    MPI_Request r;
                    int j = -1;
                    for(int k = 0; k <(int)(neighbours[dest].size()); k++){
                        if(neighbours[dest][k] == source){
                            j = k;
                            break;
                        }
                    }

                    MPI_Irecv(buffer[dest][j],3,MPI_INT, source%total_processes, source*num_stations + dest, MPI_COMM_WORLD, &r);
                    requests.push_back(r);

                }
            }
            
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
            requests.clear();
            

        }

        



        //spawn trains for this tick and place in buffer of that station
        for(int l = 0; l < 3; l++){
            char line = lines[l];
            vector<string> station_line = station_lines.at(line);
            int start = stationId[station_line.front()];
            int end = stationId[station_line.back()];

            if(Num_trains[lines[l] ] > 0){
                Num_trains[lines[l] ]-=1;
                buffer[start][6+l][0] = id;
                buffer[start][6+l][1] = l;
                buffer[start][6+l][2] = 0;  //place in start station of line
                id++;
            }

            if(Num_trains[lines[l] ] > 0){
                Num_trains[lines[l] ]-=1;
                buffer[end][6+l][0] = id;
                buffer[end][6+l][1] = l;
                buffer[end][6+l][2] = station_line.size() -1;  //place in last station of line
                id++;
            }
        }


        //sort buffers and put in queue of respective platform
        for(size_t i = mpi_rank; i < num_stations; i+=total_processes){
            insertionSort(buffer[i]);
            for(int j = 0; j < 9; j++){
                if(buffer[i][j][0] == -1)continue;
                int l = buffer[i][j][1];

                if(l < 0 || l >2)cerr<<"hi"<<endl;

                int id = buffer[i][j][0];
                int station_line_index = buffer[i][j][2];
                char line = lines[l];
                
                
                vector<string> station_line = station_lines.at(line);
                int dir;
                if(station_line_index == 0)dir = 1;
                else if(abs(station_line_index) == (int)(station_line.size()-1))dir = -1;
                else dir = station_line_index/abs(station_line_index);
                
                int next_station_line_index = abs(station_line_index) + dir;


                string next_station_name = station_line[next_station_line_index];
                int next_station_index = stationId[next_station_name];

                int temp = -1;
                for(int k = 0; k < (int)(neighbours[i].size()); k++){
                    if(neighbours[i][k] == next_station_index){
                        temp = k;
                        break;
                    }
                }
                
                queues[i][temp].push_back(make_tuple(id, l, next_station_line_index*dir));
                buffer[i][j][0] = buffer[i][j][1] = buffer[i][j][2] = -1;
            }
        }


        // update time left of train in platform or push it into link if possible
        for(size_t i = mpi_rank; i < num_stations; i+=total_processes){
            for(int j = 0; j < int(neighbours[i].size()); j++){
                if(platforms[i][j][0] == -1)continue;
                if(platforms[i][j][3] > 0){
                    platforms[i][j][3] -=1;
                } 
                if(platforms[i][j][3] == 0 && links[i][j][0] == -1){
                    links[i][j][0] = platforms[i][j][0];
                    links[i][j][1] = platforms[i][j][1];
                    links[i][j][2] = platforms[i][j][2];
                    links[i][j][3] = mat[i][neighbours[i][j]];
                    platforms[i][j][0] = platforms[i][j][1] = platforms[i][j][2] = platforms[i][j][3] = -1;  
                }
            }
        }

        // push one from queue into platfrom if possible
        for(size_t  i = mpi_rank; i < num_stations; i+=total_processes){
            for(int j = 0; j < int(neighbours[i].size()); j++){
                if(queues[i][j].size()>0 && platforms[i][j][0] == -1){
                    tuple<int,int,int> t = queues[i][j].front();
                    queues[i][j].pop_front();
                    platforms[i][j][0] = get<0>(t);
                    platforms[i][j][1] = get<1>(t);
                    platforms[i][j][2] = get<2>(t);
                    platforms[i][j][3] = gens[i][j].next(get<0>(t));
                }
            }
        }


        string ans = "";

        if(t >= ticks - num_ticks_to_print){
            //construct string for this tick 
            for(size_t i = mpi_rank; i < num_stations; i+= total_processes){
                for(int j = 0; j < (int)(neighbours[i].size()); j++){
                    
                    if(links[i][j][0] != -1){
                        string temp = "";
                        int id = links[i][j][0];
                        int l = links[i][j][1];
                        int dest = abs(links[i][j][2]);
                        int dir;
                        if(dest == 0)dir = -1;
                        else dir = dest/links[i][j][2];
                        int source = dest-dir;
                        char line = lines[l];
                        vector<string> station_line = station_lines.at(line);
                        if(source < 0 || source >= (int)(station_line.size()))cerr<<i<<" "<<source<<endl;
                        if(dest < 0 || dest >= (int)(station_line.size()))cerr<<i<<" "<<dest<<endl;

                        temp += lines[l] + to_string(id) + "-" + station_line[source] +"->" + station_line[dest]+" ";
                        ans += temp;
                    }

                    if(platforms[i][j][0] != -1){
                        int id = platforms[i][j][0];
                        int l = platforms[i][j][1];
                        string temp = "";
                        if(lines[l] != 'y' && lines[l] != 'g' && lines[l] != 'b'){
                            cerr<<l<<endl;
                        }
                        temp += lines[l]+to_string(id)+"-"+station_names[i]+"% ";
                        ans += temp;
                    }

                    for(const auto& t: queues[i][j]){
                        int id = get<0>(t);
                        int l = get<1>(t);
                        string temp = "";
                        temp += lines[l]+to_string(id)+"-"+station_names[i]+"# ";
                        ans+=temp;
                        
                    }
                }
            }
        }


        if(t >= ticks - num_ticks_to_print){

            //gather all strings from all process into rank 0 and sort it and print if needed
            int local_len = ans.size();


            if(mpi_rank == 0){
                
                for(size_t i = 1; i < total_processes; i++){
                    int str_len;
                    MPI_Recv(&str_len, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
                    vector<char> buff(str_len);
                    MPI_Recv(buff.data(), str_len, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    ans += string(buff.begin(), buff.end());
                }

                vector<string> strs;
                string temp = "";
                for(int i = 0; i < (int)(ans.size()); i++){
                    if(ans[i] == ' '){
                        strs.push_back(temp);
                        temp = "";
                    }else{
                        temp += ans[i];
                    }
                }

                sort(strs.begin(), strs.end());

                cout<<t<<": ";
                for(auto s: strs)cout<<s<<" ";
                cout<<endl;


            }else{
                MPI_Send(&local_len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Send(ans.c_str(), local_len, MPI_CHAR, 0,1,MPI_COMM_WORLD);

            }


        }
    }

}