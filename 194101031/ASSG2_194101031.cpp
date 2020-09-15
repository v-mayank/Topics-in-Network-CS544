#include <bits/stdc++.h>

using namespace std;


bool randGenProb(float prob, int seed)                // function that generat random 0 and 1 according to given prob.
{
    srand(time(0)+seed);
    float no = rand() % 100 + 1;

    if(no < prob*100)
        return 1;
    else
        return 0;
    
}

//class for INQ scheduling.
class INQswitch
{
    int numOfPorts,bufferSize, maxTime , droppedPkt,pktdpt,totalpkt; // Variables to store value of number of ports , 
    float genProb,stdDeviation;                                      // generation prob, standered deviation maxtime and etc.   
    vector< vector <float>> pktLog;             //vector to store packet logs of departing packets.
    float*** inputbuffer ;                      //input buffer for each input port.
    public:
    INQswitch(int numofport,int buffersize, int time, float pktgenprob) //class construtor to store values in variables. 
    {
        numOfPorts = numofport;             //no. of ports.
        bufferSize = buffersize;            //buffer size.
        maxTime = time;                     //max time.
        genProb = pktgenprob;               
        droppedPkt = 0;
        pktdpt = 0;
        totalpkt = 0;
        stdDeviation =0;                    //std deviation.
        inputbuffer = new float** [numOfPorts];
        for(int i=0;i<numOfPorts;i++)
        {
            inputbuffer[i]= new float* [bufferSize];
            for(int j=0;j<bufferSize;j++)
            {
                inputbuffer[i][j]= new float [4];  
                inputbuffer[i][j][0]=-1;          //arrival time.  
                inputbuffer[i][j][1]=-1;            //start time.
                inputbuffer[i][j][2]=-1;            //source port.
                inputbuffer[i][j][3]=-1;            //destination port.
            }
        }
        //inputBuffer = inputbuffer;
        
    }
    
    void departPacket(int minport,int minpos,int currTime)
    {
        //cout<<"line 5.\n";
        vector <float> log;
        //inputbuffer[minport][minpos].setDeparture(currTime);
        //cout<<"src: "<<inputbuffer[minport][minpos][2]<<"\tdest: "<<inputbuffer[minport][minpos][3]<<"\tarr: "<<inputbuffer[minport][minpos][0]<<"\tstrt: "<<inputbuffer[minport][minpos][1]<<"\tdep: "<<currTime<<endl;
        log.push_back(inputbuffer[minport][minpos][2]);         //source port
        log.push_back(inputbuffer[minport][minpos][3]);         //destination port
        log.push_back(inputbuffer[minport][minpos][0]);         //arrival time
        log.push_back(inputbuffer[minport][minpos][1]);         //start time
        log.push_back(float(currTime));                         //departure time
        log.push_back(float(currTime)-inputbuffer[minport][minpos][0]);     //delay
        //cout<<"log   "<<log[0]<<"  "<<log[1]<<"  "<<log[2]<<"  "<<log[3]<<"  "<<log[4]<<"  "<<log[5];
        pktLog.push_back(log);
        inputbuffer[minport][minpos][0]=-1;             //resetting the input buffer
        inputbuffer[minport][minpos][1]=-1;             //removing this packet from input buffer.
        inputbuffer[minport][minpos][2]=-1;
        inputbuffer[minport][minpos][3]=-1;
        pktdpt++;
    }
    void generatePacket(int currTime)
    {
        // phase 1 Packet generation.
            for(int j = 0; j<numOfPorts && currTime< maxTime;j++)
            {
                bool num = randGenProb(genProb,j);          //generate  0 or 1 .
                if(num)                                     //if 1 then packet is generated.
                {
                    //cout<<"line 1.\n";
                    totalpkt++;
                    int k ;
                    srand(time(0)+j);
                    int destinationPort = rand() % numOfPorts;  //destination port generates randomly.
                    float starttime = rand() % 10 + 1;          //start time between 0.001 and0.01
                    //cout<<"src: "<<j<<"\tdest: "<<destinationPort<<"\n";
                    for(k=0;k<bufferSize;k++){          //checking if buffer is free or not.
                        if(inputbuffer[j][k][0] == -1)  //if not packet is dropped.
                            break;
                    }
                        
                    if(k<bufferSize)
                    {
                        //cout<<"line 1.2\tk :"<<k<<endl;
                        inputbuffer[j][k][0]=currTime;              //storing values.
                        inputbuffer[j][k][1]=currTime+(starttime/1000);
                        inputbuffer[j][k][2]=j;                     //source.
                        inputbuffer[j][k][3]=destinationPort;
                    }    
                    else 
                        droppedPkt++;
                }    
                
            }
    }
    void simulation()
    {
        
        for(int currTime = 0; currTime < (maxTime+20); currTime++ )
        {
            
            generatePacket(currTime);               //packet generation phase.


            for(int outport = 0; outport<numOfPorts ; outport++)
            {
                //cout<<"line 2.\n";
                
                int minpos = -1,minport = -1;
                for(int inport = 0; inport<numOfPorts ; inport++)
                {
                    //vector <float[2]> reqForSamePort;
                    //cout<<"line 3.\n";
                    for(int index = 0;index<bufferSize;index++)
                    {
                        float minTime = 999999;
                        if(outport==inputbuffer[inport][index][3])
                        {
                            //reqForSamePort.push_back(inputbuffer[inport][index][1]);
                            if(minTime > inputbuffer[inport][index][1])         //finding min start time in buffer of a port.
                            {
                                minTime = inputbuffer[inport][index][1];
                                minpos = index;                     //storing its position and port no.
                                minport = inport;
                                //cout<<"line 4.\n";
                            }    
                        }
                    }
                    
                }
                //cout<<"minport: "<<minport<<"\tminpos: "<<minpos<<endl;
                if(minport!=-1&&minpos!=-1)
                {
                    //cout<<"minport: "<<minport<<"\tminpos: "<<minpos<<"\tdest: "<<inputbuffer[minport][minpos][3]<<endl;
                    departPacket(minport,minpos,currTime);      //depart that packet.
                }    
            }    

        }
    }
    void display()
    {
        simulation();
        //cout<<"src port\tdest port\tarrival\tstart\tdep\tdelay\n";
        float sum = 0;
        for(int row=0;row<pktLog.size();row++)
        {
            ////cout<<"line 10.\n";
            sum = sum + pktLog[row][5];

          /*  for(int col=0;col<pktLog[row].size();col++)
            {
                
                cout<<pktLog[row][col]<<"   ";
            }
            cout<<endl;*/
        }
        double avgdelay = sum/pktdpt;
        float linkUtilization = (pktdpt)/(pktLog.back()[4] * numOfPorts * 1.0);
        for(int itr=0;itr<pktLog.size();itr++)
        {
            stdDeviation = stdDeviation + ((pktLog[itr][5] - avgdelay)*(pktLog[itr][5] - avgdelay));  
        }
        stdDeviation = sqrt(stdDeviation/pktLog.size());
        char queueType[5]="INQ";
        ofstream fout;
        fout.open("out.txt",ios::app);
        //fout<<"Number of ports\t Probability\t QueueType\t Average Delay\t Standard Deviation of Delay\t Link Utilisation";
        fout << numOfPorts << "\t" << genProb << "\t" << queueType << "\t\t" << avgdelay << "\t\t" << stdDeviation << "\t\t" << linkUtilization<<endl;
        fout.close();
        cout<<"avg delay: "<<sum/pktdpt<<"\tpkct depart: "<<pktdpt<<"\t pckt dropped: "<<droppedPkt<<"\ttotal pkts: "<<totalpkt<<endl;

        cout<<"\n\nN = "<<numOfPorts<<"\tprob = "<<genProb<<"\tINQ\tAvg delay = "<<avgdelay<<"\tStd dev = "<<stdDeviation<<"\tAvg link Util."<<linkUtilization<<endl;
    }
};
//vector to store output of KOUQ scheduling.
vector<vector<float>> output;
//function to create output
void create_output(int N)             
{
    for(int i=0;i<N;i++)
    {
        vector<float> v;
        output.push_back(v);
    }
}
//function to packet generation for KOUQ scheduling.
double packet_generation(int N,int B,float prob,int max_time,int K,int t)
{
   
    vector<int> packet_at_dest;
    double packet_drop=0;
        packet_at_dest.assign(N,0);
        for(int i=0;i<N;i++)
        {
           float ran_prob = ((rand()+rand())*rand())%100;
          // double ran_prob=((double) rand() / (RAND_MAX));
            ran_prob /= 100;
            if(!(ran_prob>prob))
            {
                int dest=(rand()%N);
                packet_at_dest[dest]++;
                if(!(packet_at_dest[dest]<K) && !(output[dest].size()<B))
                packet_drop++;
                else
                {
                    output[dest].push_back(t+ran_prob);   
                }
                
            }

        }
        return packet_drop;
}
//function to KOUQ scheduling.
void KOUQ_scheduling(int N,int B,float prob,int max_time)
{
    create_output(N);
    int K=N*1.0;
    double packet_drop;
    double link_util=0;
    vector<double> delay;
    for(int t=0;t<max_time;t++)
    {
        packet_drop=packet_generation(N,B,prob,max_time,K,t);
        for(int i=0;i<N;i++)
        {
             if(output[i].size()!=0 )
            {
                delay.push_back(t-output[i].front()+1);
                output[i].erase(output[i].begin());
                link_util++;
            }
        }
    }
    double avg_delay=0;
    for(int i=0;i<delay.size();i++)
        avg_delay+=delay[i];
    avg_delay=avg_delay/link_util;
    float stand_dev = 0;
    for (int i = 0; i < delay.size(); i++) {
        stand_dev += (delay[i] - avg_delay) * (delay[i] - avg_delay); 
    }
    stand_dev = sqrt(stand_dev / delay.size());
    double link=link_util/(N*max_time*1.0);
    double drop=packet_drop/(N*max_time*1.0);
    cout<<"KOUQ:"<<endl;
    cout<<"DELAY:"<<avg_delay<<endl;
    cout<<"LINK UTILISATION:"<<link<<endl;
    cout<<"PACKET DROP PROBABILITY:"<<drop<<endl;
    cout<<"Standard deviation :"<<stand_dev<<endl;
    ofstream fout;
    fout.open("out.txt",ios::app);
        //fout<<"Number of ports\t Probability\t QueueType\t Average Delay\t Standard Deviation of Delay\t Link Utilisation";
        fout << N << "\t" << prob << "\t" << "KOUQ" << "\t\t" << avg_delay << "\t\t" << stand_dev << "\t\t" << link<<endl;
        fout.close();

    cout<<"\nN = "<<N<<"\tprob = "<<prob<<"\tKOUQ\tAvg delay = "<<avg_delay<<"\tStd dev = "<<stand_dev<<"\tAvg link Util."<<link<<endl;

}


//structure to store information of packet
typedef struct packet {
    bool is_generated; //if packet is generated or not
    int ip; //input port
    int op; //output port
    double gen_time; //geneartion tym of packet
    double comp_time; //completion time of transmission
    //We will calculate delay by differenece of above two
}
packet;

//Function to Generate accept Packet
packet generate_packet(int ip, int t, float p, int numOfports) {
    packet pkt;
    pkt.is_generated = randGenProb(p,t); //genearte accept random no. between 0 and 1
    int op = rand() % numOfports; //selecting accept random output port 
    double time = t + (double)((rand() % 10 + 1) / 1000.0); //generation time selected randomly between t+0.001 to t+0.01
    //setting claculated params to packet parameteers
    pkt.ip = ip;
    pkt.gen_time = time;
    pkt.op = op;
    //returning packet
    return pkt;
}
//Function to islip scheduling.
void Islip(int numOfports, float pktGenProb, int maxTime ) {
    //declaration
    char queueType[20] = "islip";
    int delay = 0; //delay variable that stores delay sum over all simulation time
    int pktTransmitted = 0; //to store number of transmitted packet overall
    int genPkt = 0; //to store number of genearetd packet
    vector < int > delay_arr; //delay vector to claculate std. deviation
    queue < packet > ipQueue[numOfports][numOfports];
    int acceptPointer[numOfports] = { //Accept pointer array
        0
    };
    int grantPointer[numOfports] = { //Grant Pointer array
        0
    };
    srand(time(0));
    for (int t = 0; t < maxTime; t++) {
        //cout<<"Round : "<<t<<endl;
        int request[numOfports][numOfports]; //Array to maintain requests from input-port
        int grant[numOfports] = {-1}; //Array to maintain Grant-requests
        int accept[numOfports] = {-1}; //Array to Accept Requests
        bool isOPconnected[numOfports] = {false}; //Array to check whether the give output-port connected or not
        bool isIPconnected[numOfports] = {false}; //Array to check whether given input port has established connection
        //generate packets for each input-port in each time-slot
        for (int ip = 0; ip < numOfports; ip++) 
        {
            packet p = generate_packet(ip, t, pktGenProb, numOfports); //generate packet
            if (p.is_generated) 
            {   //check if packet is generated to insert packet in queue
                genPkt++; //Maintain count of generated packets
                ipQueue[ip][p.op].push(p); //Push the packet in specified queue
            }
        }
        int conn = 1, iter = 0; 
        while (conn > 0)
        {
            conn = 0;
            //Build Requests
            // cout<<"Iter : "<<iter<<endl;
            // cout<<"Requested Connections : "<<endl;
            for (int ip = 0; ip < numOfports; ip++) 
            {
                for (int op = 0; op < numOfports; op++) 
                {
                    if (!ipQueue[ip][op].empty()&& ( !isIPconnected[ip] && !isOPconnected[op] )) 
                    {   
                        // Check whether queue has packet
                        request[ip][op] = 1; //Make request
                        //cout<<ip<<" : "<<op<<endl;
                    } 
                    else 
                        request[ip][op] = -1; 
                }
            }

            //Grant Requests
            //cout<<"Granted connections : "<<endl;
            for (int op = 0; op < numOfports; op++) 
            {
                if (!isOPconnected[op])        //Check if output-port established connection 
                {
                    int i = 0, ip = grantPointer[op];
                    for (; i < numOfports; i++) 
                    {   //Check for next avilable request
                        if (request[ip][op] == 1) break; // Break if request found
                        ip = (ip + 1) % numOfports;
                    }
                    if (i < numOfports)
                    {
                        grant[op] = ip;
                        //cout<<op<<" : "<<ip<<endl; //Grant to first avilable request
                    } 
                    else 
                        grant[op] = -1;
                }           
                
                
            }

            //Accept-phase
            //cout<<"Accepted connections : "<<endl;
            for (int ip = 0; ip < numOfports; ip++) 
            {
                if (!isIPconnected[ip]) //Check if input-port is busy
                {
                    int i = 0, op = acceptPointer[ip];
                    for (; i < numOfports; i++) 
                    { //Check for first port that granted request
                        if (grant[op] == ip) break; //Break if Grant found
                        op = (op + 1) % numOfports;
                    }
                    if (i < numOfports) 
                    {
                        if( iter == 0 )
                        {
                            acceptPointer[ip] = (op+1)%numOfports;
                            grantPointer[op] = (ip+1)%numOfports;
                        }
                        
                        //cout<<ip<<" : "<<op<<endl;
                        isOPconnected[op] = true; // Mark outport as established connection
                        isIPconnected[ip] = true; // Mark inputport as establishe connection
                        packet p = ipQueue[ip][op].front(); //store the front packet in the queue
                        p.comp_time = t + 1; //store completion time
                        delay += (int) p.comp_time - (int) p.gen_time; //calculate total delay till now
                        delay_arr.push_back((int) p.comp_time - (int) p.gen_time); //maintain delay of each packet
                        ipQueue[ip][op].pop(); //remove packet from queue
                        conn++; //maintain no.of connections established in each iteration
                        pktTransmitted++; //Count no.of packets transmitted
                    }
                }
                
            }
            iter++;
        } 
    }
    float AvgDelay = delay / (pktTransmitted * 1.0); //Calculate average delay
    // Calculate Standard-Deviation
    float stand_dev = 0;
    for (int i = 0; i < delay_arr.size(); i++) {
        stand_dev += (delay_arr[i] - AvgDelay) * (delay_arr[i] - AvgDelay); 
    }
    stand_dev = sqrt(stand_dev / delay_arr.size());
    // calculate link utilisation
    float linkUtill = pktTransmitted / (numOfports * maxTime * 1.0);
    ofstream fout;
    fout.open("out.txt", ios::app);
    fout << numOfports << "\t" << pktGenProb << "\t" << queueType << "\t\t" << AvgDelay << "\t\t" << stand_dev << "\t\t" << linkUtill<<endl;
    fout.close();

    cout<<"\n\nN = "<<numOfports<<"\tprob = "<<pktGenProb<<"\tISLIP\tAvg delay = "<<AvgDelay<<"\tStd dev = "<<stand_dev<<"\tAvg link Util."<<linkUtill<<endl;
    // fout.open("delay_onprob.txt",ios::app);
    // fout<<numOfports<<"\t"<<AvgDelay<<endl;    
    // fout.close();
    // fout.open("link_onprob.txt",ios::app);
    // fout<<numOfports<<"\t"<<linkUtill<<endl;    
    // fout.close();
}



int main(int argc, char* argv[])
{
    int portcount = 8;  //number of ports with default value.
    int buffersize = 4; //buffer size.
    float pktGenProb = 0.5 ;    //probability of packet generation.
    const char* queueType = "INQ";  //queue scheduling type.
    int maxTimeSlot = 10000;    //maximum time slots.
    float knockout = 0.6;       //knockout out values.
    //cout<<argc<<endl;
    for(int i = 1 ; i < argc ; i++)             //parsing the arguments passed by user.
    {
        //cout<<argv[i]<<endl;
        char* temp = argv[i] ;
        char* token;
        token = strtok(temp,"_");               //breaking arguments into tokens before and after "_"(underscore).
        //cout<<token<<endl;
        if(strcmp(token,"N")==0)
        {
            token = strtok(NULL,"_");           //though "_" is look like smiley(emoticon) but it is a delimiter.
            //cout<<token<<endl;
            portcount = atoi(token);
            //cout<<portcount<<endl;
        }
        if(strcmp(token,"B")==0)                //B_size 
        {
            token = strtok(NULL,"_");
            //cout<<token<<endl;
            buffersize = atoi(token);
            //cout<<buffersize<<endl;
        }
        if(strcmp(token,"P")==0)                //P_prob
        {
            token = strtok(NULL,"_");
            //cout<<token<<endl;
            pktGenProb = atof(token);
            if(pktGenProb<0||pktGenProb>1)
            {
                cout<<"Enter valid options.\n";
                return -1;
            }
            
        }
        if(strcmp(token,"Q")==0)                //Q_queuetype
        {
            token = strtok(NULL,"_");
            //cout<<token<<endl;
            queueType = token;
            //cout<<queueType<<endl;
        }
        if(strcmp(token,"K")==0)                //k_knocoutvalue
        {
            token = strtok(NULL,"_");
            //cout<<token<<endl;
            knockout = atof(token);
            if(knockout<0||knockout>1)
            {
                cout<<"Enter valid options.\n";
                return -1;
            }
            
            //cout<<knockout<<endl;
        }
        if(strcmp(token,"T")==0)                //T_maxtime
        {
            token = strtok(NULL,"_");
            //cout<<token<<endl;
            maxTimeSlot = atoi(token);
            //cout<<maxTimeSlot<<endl;
        }        

    }
    cout<<portcount<<endl;
    cout<<buffersize<<endl;
    cout<<maxTimeSlot<<endl;
    cout<<pktGenProb<<endl;
   
   if(strcmp(queueType,"INQ")==0)       //checking type of queue to run scheduling. INQ is default.
   {
    INQswitch s1(portcount,buffersize,maxTimeSlot,pktGenProb);
    s1.display();
   }
   else if(strcmp(queueType,"KOUQ")==0)
   {
    KOUQ_scheduling(portcount,buffersize,pktGenProb,maxTimeSlot); 
   }
   else if(strcmp(queueType,"ISLIP")==0)
   {
       Islip(portcount,pktGenProb,maxTimeSlot);
   }
   else
   {
       
       cout<<"Bhai shi naam daal na\n";
   }
   
    return 0;
}