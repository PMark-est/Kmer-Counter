/**
 * MEMORY LEAK FREE PROGRAM
 * In case you're a poor soul debugging or analysing this mess I'll give a very brief overview of how this works
 * First of all the whitespace is removed from the fasta file that is being worked on
 * Secondly kmers are read as characters into a string that is converted to binary
 * Thirdly the kmer is stored as a unsigned long in a chained hashtable. Theres no hash function since each kmer
 * has its own respective binary representation(a = 00, c = 01, g = 10, t = 11). I.E acg = 00 01 10 = 6; this is stored
 * in slot 6 of the hashtable. You may say that cg, acg, aacg, aaacg, etc are the same. This is true
 * but the aforementioned scenario is impossible since there is a K value that the user inputs which determines how long
 * the K-mer is.
 */

#include <iostream>
#include <thread>
#include <filesystem>
#include <vector>
#include <bitset>
#include <fstream>
#include <cmath>

#define MAX_SIZE 64 // Biggest kmer in bits, 32-mer
#define MB 1048576.0
size_t kmerMax = 0; // number of kmer combinations(4^k)
std::string deleteMetadata = "sed -i \'/>/d\' "; // sed -i '/>/d'
std::string deleteEmptyLines = "sed -i \'/^$/d\' "; // sed -i '/^$/d'

size_t hash_c_string(const char* p, size_t size) {
    size_t result = 0;
    const size_t prime = 31;
    for (size_t i = 0; i < size; ++i) {
        result = p[i] + (result * prime);
    }
    return result;
}

char* copy(const char *string, size_t size){
    char *copy = new char[size+1];
    for (int i = 0; i < size; i++) {
        copy[i] = string[i];
    }
    copy[size] = '\0';
    return copy;
}

bool equals(const char *string1, const char *string2, size_t size){
    if (string1 == nullptr || string2 == nullptr) return false;
    for (int i = 0; i < size; i++){
        if (string1[i] != string2[i]) return false;
    }
    return true;
}

class HashTable{
public:
    char **keys;
    char *vals;
    int length;
    explicit HashTable(int size){
        keys = new char*[size];
        for (int i = 0; i < size; ++i) {
            keys[i] = nullptr;
        }
        vals = new char[size];
        length = size;
    }
    ~HashTable(){
        for (int i = 0; i < length; ++i) {
            delete[] keys[i];
        }
        delete[] keys;
        delete[] vals;
    }
    void put(char *key, size_t size, char val) const{
        size_t hash = hash_c_string(key, size) % length;
        int threshold = 0;
        while (keys[hash] != nullptr){
            if (threshold == 10) return;
            if(equals(keys[hash], key, size)) return;
            hash++;
            hash %= length;
            threshold++;
        }
        keys[hash] = copy(key, size);
        vals[hash] = val;
    }
    char get(char *string, size_t size) const{
        size_t hash = hash_c_string(string, size) % length;
        int threshold = 0;
        while (!equals(keys[hash], string, size)){
            if (threshold == 10) return '\0';
            hash++;
            hash %= size;
            threshold++;
        }
        return vals[hash];
    }
};

class Node{
public:
    size_t resOccurences{};
    size_t susOccurences{};
    size_t fileNr{};
    size_t data{};
    Node *next{};
    Node()= default;
    Node(Node &other): resOccurences(other.resOccurences), susOccurences(other.susOccurences), fileNr(other.fileNr), data(other.data){
        next = nullptr;
    }
    explicit Node(unsigned long data){
        this->data = data;
        this->next = nullptr;
    }
};

class LinkedList{
public:
    Node *head{};
    LinkedList()= default;
    ~LinkedList(){
        Node *tmp = head;
        while(head){
            head = head->next;
            delete tmp;
            tmp = head;
        }
    }
    void push(Node *node){
        // Linked list is empty
        if (head == nullptr){
            head = node;
            node->next = nullptr;
            return;
        }
        // Prepends new node to linked list
        Node *tmp = head;
        head = node;
        head->next = tmp;
    }

    void push(unsigned long data, size_t fileNr, bool isRes, size_t *nodeCount){
        // Linked list is empty
        if (head == nullptr){
            Node *newNode = new Node(data);
            head = newNode;
            newNode->fileNr=fileNr;
            if (isRes) newNode->resOccurences++;
            else newNode->susOccurences++;
            (*nodeCount)++;
            return;
        }
        Node *tmp = head;
        // Linked list's first element is equal to given data
        if (tmp->data == data){
            if (tmp->fileNr != fileNr){
                if (isRes) tmp->resOccurences++;
                else tmp->susOccurences++;
            }
            tmp->fileNr = fileNr;
            return;
        }
        // Traverse linked list until it finds a node with the given data or until it reaches the end of the list.
        while (tmp->next != nullptr){
            if (tmp->data == data){
                if (tmp->fileNr != fileNr){
                    if (isRes) tmp->resOccurences++;
                    else tmp->susOccurences++;
                    tmp->fileNr = fileNr;
                }
                return;
            }
            tmp = tmp->next;
        }
        // Prepends new node to linked list
        tmp = head;
        head = new Node(data);
        head->next = tmp;
        if (isRes) head->resOccurences++;
        else head->susOccurences++;
        head->fileNr = fileNr;
        (*nodeCount)++;
    }

    void print() const{
        Node *tmp = head;
        while (tmp != nullptr){
            std::cout << tmp->data << "(Res: " << tmp->resOccurences << ", Sus: " << tmp->susOccurences << ")->";
            tmp = tmp->next;
        }
        std::cout << "null\n";
    }
};

size_t convertToDecimal(char **binary, size_t k){
    size_t decimalNumber = 0;
    auto pow = (size_t) std::pow(2, k-1);
    for (int i = 0; i < k; i++) {
        decimalNumber += ((*binary)[i] - 48) * pow;
        pow /= 2;
    }
    return decimalNumber;
}

// Copies elements from smaller array to bigger array
void resizeBigger(LinkedList *&nodes, size_t &listSize){
    size_t newSize = listSize << 1;
    auto *newArr = new LinkedList[newSize];
    int i = 0;
    Node *tmp;
    for ( ; i < listSize; i++){
        if (nodes[i].head) {
            tmp = nodes[i].head;
            while(tmp){
                size_t index = tmp->data % newSize;
                newArr[index].push(new Node(*tmp));
                tmp = tmp->next;
            }
        }
    }
    delete[] nodes;
    nodes = newArr;
    listSize = newSize;
}

void writeToFile(LinkedList **lists, size_t *listSizes, const size_t threadCount, const size_t k, int fileAmount){
    std::ofstream kmersFile;
    kmersFile.open("counts.csv");
    kmersFile << k << "-mer(convert to binary (2*k) to get nucleotides; 00=A,01=C,10=G,11=T),res,sus\n";
    auto *sortingTable = new LinkedList[fileAmount];
    std::cout << "started writing\n";
    if (threadCount == 1){
        LinkedList *list = *lists;
        Node *tmp;
        Node *tmpNext;
        for (int i = 0; i < listSizes[0]; i++){
            if(list[i].head){
                // copies address to new table and sets old address location to null
                // this makes the destructor not go apeshit
                tmp = list[i].head;
                while(tmp){
                    list[i].head = nullptr;
                    tmpNext = tmp->next;
                    list[i].head = tmpNext;
                    sortingTable[tmp->resOccurences + tmp->susOccurences - 1].push(tmp); // res + sus is always at least 1
                    tmp = tmpNext;
                }
            }
        }
        int total = 0;
        for (int i = fileAmount-1; i > -1; i--) {
            Node *node = sortingTable[i].head;
            int p = 0;
            while (node){
                kmersFile << node->data << ", " << node->resOccurences << "," << node->susOccurences << "\n";
                node = node->next;
                p++;
            }
            total += p;
        }
   }
    else{
        size_t tempSize;
        LinkedList *tempListPtr;
        // sort lists by size in ascending order
        for (int i = 0; i < threadCount; i++) {
            for (int j = 0; j < threadCount; j++) {
                if(listSizes[i] < listSizes[j]){
                    tempSize = listSizes[i];
                    tempListPtr = lists[i];
                    listSizes[i] = listSizes[j];
                    listSizes[j] = tempSize;
                    lists[i] = lists[j];
                    lists[j] = tempListPtr;
                }
            }
        }
        Node *tmp;
        Node *tmp2;
        Node *prev;
        size_t index;
        for (int i = 0; i < threadCount; i++) {
            for (int j = 0; j < listSizes[i]; j++) {
                tmp = lists[i][j].head;
                while(tmp){
                    for (int l = i+1; l < threadCount; l++) {
                        index = tmp->data % listSizes[l];
                        prev = lists[l][index].head;
                        if (!prev) continue;
                        else{
                            if (prev->data != tmp->data) {
                                tmp2 = prev->next;
                            }
                            else{
                                tmp->susOccurences += prev->susOccurences;
                                tmp->resOccurences += prev->resOccurences;
                                delete lists[l][index].head;
                                lists[l][index].head = prev->next;
                                continue;
                            }
                        }
                        while(tmp2 && tmp2->data != tmp->data){
                            tmp2 = tmp2->next;
                            prev = prev->next;
                        }
                        if (tmp2){
                            tmp->susOccurences += tmp2->susOccurences;
                            tmp->resOccurences += tmp2->resOccurences;
                            prev->next = tmp2->next;
                            delete tmp2;
                        }
                    }
                    tmp2 = tmp->next;
                    sortingTable[tmp->resOccurences+tmp->susOccurences-1].push(tmp);
                    tmp = tmp2;
                    //delete lists[i][j].head;
                    lists[i][j].head = nullptr;
                    lists[i][j].head = tmp;
                }
            }
        }
        for (int i = fileAmount-1; i > -1; i--) {
            tmp = sortingTable[i].head;
            while(tmp){
                kmersFile << tmp->data << ", " << tmp->resOccurences << "," << tmp->susOccurences << "\n";
                tmp = tmp->next;
            }
        }
    }
    std::cout << "writing done\n";
    delete[] sortingTable;
    kmersFile.close();
}

// The fileSize may be bigger than the amount of kmers since the file has newlines (\n)
void readFile(const size_t k, const size_t fileSize, const size_t *listSize, size_t *nodeCount, size_t occupiedFile, bool isRes, char **kmer, char *nucleotides, LinkedList *list){
    size_t currentNucleotideIndex = 0;
    for (int i = 0; i < 2*k; i+=2, currentNucleotideIndex++){
        switch (*nucleotides) {
            case 'a': (*kmer)[i] = 48; (*kmer)[i+1] = 48; break;
            case 'c': (*kmer)[i] = 48; (*kmer)[i+1] = 49; break;
            case 'g': (*kmer)[i] = 49; (*kmer)[i+1] = 48; break;
            case 't': (*kmer)[i] = 49; (*kmer)[i+1] = 49; break;
            case 'A': (*kmer)[i] = 48; (*kmer)[i+1] = 48; break;
            case 'C': (*kmer)[i] = 48; (*kmer)[i+1] = 49; break;
            case 'G': (*kmer)[i] = 49; (*kmer)[i+1] = 48; break;
            case 'T': (*kmer)[i] = 49; (*kmer)[i+1] = 49; break;
        }
        nucleotides++;
    }
    unsigned long val = convertToDecimal(kmer, k);
    unsigned int key = val % *listSize;
    list[key].push(val, occupiedFile, isRes, nodeCount);
    for (; currentNucleotideIndex < fileSize; currentNucleotideIndex++){
        if (*nucleotides != '\n'){
            val <<= 2;
            val %= kmerMax;
            // case 'a' is not needed as += 0 is the same as doing nothing
            switch (*nucleotides) {
                case 'c': val += 1; break;
                case 'g': val += 2; break;
                case 't': val += 3; break;
                case 'C': val += 1; break;
                case 'G': val += 2; break;
                case 'T': val += 3; break;
            }
            key = val % *listSize;
            list[key].push(val, occupiedFile, isRes, nodeCount);
        }
        nucleotides++;
    }
}

void readFiles(const size_t k, size_t *listSize, size_t *nodeCount, const size_t initialBufferSize, const std::vector<std::filesystem::directory_entry> &files, const std::vector<bool> &resistances, LinkedList **list){
    auto start = std::chrono::high_resolution_clock::now();
    char *kmer = new char[2*k];
    size_t bufferSize = initialBufferSize;
    char *buffer = new char[bufferSize];
    size_t fileNr = 1;
    for (int i = 0; i < files.size(); ++i) {
        std::string fileName = files[i].path().string();
        bool isResistant = resistances[i];
        // Remove useless lines from fna file. like lines that begin with >>> or empty lines
        if (system((deleteMetadata + fileName).c_str()) || system((deleteEmptyLines + fileName).c_str()))
            std::cout << "Removed lines from " << fileName << "\n";
        std::ifstream genomeFile(fileName);

        // Read how many bytes the file is
        genomeFile.seekg(0, std::ios::end);
        auto fileSize = genomeFile.tellg();
        if (fileSize > bufferSize){
            bufferSize = (int) fileSize << 1;
            delete[] buffer;
            buffer = new char[bufferSize];
        }
        genomeFile.seekg(0, std::ios::beg);
        if (( (int) fileSize + *nodeCount) > ((*listSize) >> 1)) resizeBigger(*list, *listSize);
        genomeFile.read(buffer, fileSize);
        readFile(k, fileSize, listSize, nodeCount, fileNr, isResistant, &kmer, buffer, *list);
        fileNr++;
    }
    // Calculate and display how long the file was read for
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << duration.count() / 1000000.0 << " ms\n";
    delete[] kmer;
    delete[] buffer;
}

HashTable* readMetadataToTable(const std::string path){
    std::ifstream metadata(path);
    int fileSize = 0;
    std::string line;
    while(std::getline(metadata, line)) {
        fileSize++;
    }
    fileSize <<= 1;
    auto resistanceTable = new HashTable(fileSize);
    metadata.clear();
    metadata.seekg(0);
    std::getline(metadata, line);
    int comma;
    char *genomeID = new char[64];
    int phenotype;
    char val;
    int i;
    while(std::getline(metadata, line)){
        std::stringstream ss(line);
        i = 0; // keeps track of the index at which a character gets placed in the array when making the genome id
        phenotype = 0; // Keeps track of how long the word for the phenotype is(resistant or susceptible)
        // because of the string stream the quotation marks are kept, so they're "resistant" or "susceptible"
        comma = 0; // Keeps track of how many commas have been seen
        while (ss >> val) {
            // gets phenotype
            if(comma == 4){
                phenotype++;
                if (phenotype == 2 && val != 'S' && val != 'R') break;
            }
            // gets genome id
            if (comma == 1){
                if(val == '"') continue;
                genomeID[i] = val;
                i++;
            }
            if (val == ',') {
                comma++;
                if (comma == 2){
                    genomeID[i-1] = '\0';
                }
                else if (comma == 5) {
                    if (phenotype == 14) resistanceTable->put(genomeID, i, 's');
                    else if(phenotype == 12) resistanceTable->put(genomeID, i, 'r');
                    break;
                }
            }
        }
    }
    delete[] genomeID;
    metadata.close();
    return resistanceTable;
}

std::string getPath(){
    std::ifstream settingsFileI("settings.txt");
    std::string path;
    if (settingsFileI){
        getline(settingsFileI, path);
    } else{
        std::cout << "no settings file found\nPlease enter the absolute path containing the fna files: ";
        std::cin >> path; // /home/marko/Git/Bio/files/
        std::ofstream settingsFileO("settings.txt");
        settingsFileO << path;
        settingsFileO.close();
    }
    settingsFileI.close();
    return path;
}



int main(int argc, char* argv[]){
    std::string path = getPath();
    unsigned int k;
    int threadCount;
    std::string bacterium;
    const auto processor_count = std::thread::hardware_concurrency();
    if (argc > 1){
        bacterium = argv[1];
        try{
            threadCount = std::stoi(argv[2]);
            std::string input;
            if (threadCount > processor_count){
                std::cout << "Can't select more than " << processor_count << " threads\n";
                return 0;
            }
            else if (threadCount > 0.5 * processor_count){
                std::cout << "****WARNING****\nYOU HAVE SELECTED " << threadCount << " THREADS TO MULTI-THREAD WITH\nARE YOU SURE YOU WANT TO CONTINUE(y/n)?: ";
                std::cin >> input;
                if (input != "Y" || input != "y"){
                    std::cout << "\nABORTED\n";
                    return 0;
                }
            }
            k = std::stoi(argv[3]);
        } catch (const std::exception& e){
            std::cout << "didn't enter a number\n";
            return 0;
        }
    } else{
        try{
            std::string input;
            std::cout << "How many threads do you want to use? You have " << processor_count << " threads to use: ";
            std::cin >> input;
            threadCount = std::stoi(input);
            if (threadCount > processor_count){
                std::cout << "Can't select more than " << processor_count << " threads\n";
                return 0;
            }
            else if (threadCount > 0.5 * processor_count){
                std::cout << threadCount << ", " << 0.5*threadCount << "\n";
                std::cout << "****WARNING****\nYOU HAVE SELECTED " << threadCount << " THREADS TO MULTI-THREAD WITH\nARE YOU SURE YOU WANT TO CONTINUE(y/n)?: ";
                std::cin >> input;
                if (input != "Y" || input != "y"){
                    std::cout << "\nABORTED\n";
                    return 0;
                }
            }
            std::cout << "What k value do you want(less than 32)? ";
            std::cin >> input;
            k = std::stoi(input);
            if (k > 31){
                std::cout << "k isn't less than 32!, K: " << k << "\n";
            }
        } catch (const std::exception& e){
            std::cout << "didn't enter a number\n";
            return 0;
        }
        while(true){
            std::cout << "-----\nTo view your folders, type \'h\'\n";
            std::cout << "What bacterium would you like to analyze? ";
            std::cin >> bacterium;
            if(bacterium=="h" || bacterium=="H") {
                for (const auto &entry: std::filesystem::directory_iterator(std::filesystem::path(path))) {
                    std::cout << entry.path().filename().string() << "\n";
                }
            } else
                break;
        }
    }
    kmerMax = static_cast<size_t>(std::pow(4, k));
    auto **lists = new LinkedList *[threadCount]; // array of linked hash tables
    auto *listSizes = new size_t[threadCount]; // sizes of the hash tables
    auto *nodeCounts = new size_t[threadCount]; // stores how many elements are in a hash table
    auto *threads = new std::thread[threadCount];
    auto *files = new std::vector<std::filesystem::directory_entry>[threadCount]; // array of vectors that hold the file names for multi-threading
    auto *resistances = new std::vector<bool>[threadCount]; // array of vectors that hold whether the file is resistant or susceptible
    int i = 0;
    int fileCount = 0;
    char *fileName;
    std::string fileNameS;
    path += bacterium;
    auto folder = std::filesystem::directory_iterator(path);
    auto table = readMetadataToTable(folder->path().filename().string());
    folder++;
    while (folder){
        i %= threadCount;

        fileNameS = folder->path().filename().string();
        fileName = copy(fileNameS.c_str(), fileNameS.length()-4);
        char resistance = table->get(fileName, fileNameS.length()-3);
        delete fileName;
        if(resistance == 'r') resistances[i].push_back(true);
        else if(resistance == 's') resistances[i].push_back(false);
        else continue;

        files[i].push_back(entry);
        i++;
        fileCount++;
        folder++;
    }
    delete table;
    size_t fileSize;
    for (i = 0; i < threadCount; i++) {
        std::ifstream genomeFile(files[i].back().path().string());
        genomeFile.seekg(0, std::ios::end);
        fileSize = genomeFile.tellg();
        fileSize <<= 1;
        genomeFile.close();
        listSizes[i] = fileSize;
        lists[i] = new LinkedList[fileSize];
        nodeCounts[i] = 0;

        threads[i] = std::thread(readFiles, k, &listSizes[i], &nodeCounts[i], fileSize, files[i], resistances[i], &lists[i]);
    }
    for (i = 0; i < threadCount; i++)
        threads[i].join();

    writeToFile(lists, listSizes, threadCount, k, fileCount);
    std::cout << "Finished\n";
    // free memory
    for (int j = 0; j < threadCount; j++) {
        delete[] lists[j];
    }
    delete[] lists;
    delete[] listSizes;
    delete[] nodeCounts;
    delete[] threads;
    delete[] files;
    delete[] resistances;
    return 0;

}
