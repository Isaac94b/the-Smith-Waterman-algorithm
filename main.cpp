#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "termcolor\termcolor.hpp"

using namespace std;

// Struct for sequence alignment properties
typedef struct SM_Prop {
    string First_Seq;
    string Second_Seq;
    int S_Matrix[50][50];
} SM_Prop;

// Struct for matrix rules (match, mismatch, gap)
typedef struct Matrix_Rules {
    int Match;
    int Mismatch;
    int Gap;
} M_Rules;

// Function prototypes
int Max(int HV, int VV, int DV);
void Initialization();
void Scoring();
void traceBackANDPrint();

// Global variables
M_Rules MMG;
SM_Prop FSS;
vector<int> PathT;
vector<int> PathT2;
vector<string> Seq;
vector<string> Seq2;
vector<string> Seqq;
vector<string> Seqq2;

int main() {
    // Input match, mismatch, and gap values
    cout << termcolor::yellow << "\n Please enter the Match value : " << termcolor::reset;
    cin >> MMG.Match;

    cout << termcolor::yellow << "\n Please enter the Mismatch value : " << termcolor::reset;
    cin >> MMG.Mismatch;

    cout << termcolor::yellow << "\n Please enter the Gap value : " << termcolor::reset;
    cin >> MMG.Gap;

    // Input sequences (commented out to use default sequences for testing)
    // cout << " Please enter the first sequence : ";
    // cin >> FSS.First_Seq;
    // cout << " Please enter the second sequence : ";
    // cin >> FSS.Second_Seq;

    // Default sequences for testing
    FSS.First_Seq = "ACAGTC";
    FSS.Second_Seq = "CATTGC";

    // Alignment process
    Initialization();
    Scoring();
    traceBackANDPrint();

    // Wait for user input before exiting
    getchar();
    cout << "\n";
    getchar();
    return 0;
}

// Function to calculate the maximum of three values
int Max(int HV, int VV, int DV) {
    vector<int> Tmpr = {HV, VV, DV};
    sort(Tmpr.begin(), Tmpr.begin() + Tmpr.size());
    return Tmpr[2];
}

// Function for initialization of the scoring matrix
void Initialization() {
    for (int i = 0; i <= FSS.First_Seq.size(); i++) {
        for (int j = 0; j <= FSS.Second_Seq.size(); j++) {
            if (i == 0 && j == 0) {
                FSS.S_Matrix[i][j] = 0;
            } else if (i == 0 && j > 0) {
                FSS.S_Matrix[0][j] = FSS.S_Matrix[0][j - 1] + MMG.Gap;
            } else if (i > 0 && j == 0) {
                FSS.S_Matrix[i][0] = FSS.S_Matrix[i - 1][0] + MMG.Gap;
            } else {
                if (FSS.First_Seq[i - 1] == FSS.Second_Seq[j - 1]) {
                    FSS.S_Matrix[i][j] = Max(FSS.S_Matrix[i - 1][j] + MMG.Gap, FSS.S_Matrix[i][j - 1] + MMG.Gap,
                                             FSS.S_Matrix[i - 1][j - 1] + MMG.Match);
                } else {
                    FSS.S_Matrix[i][j] = Max(FSS.S_Matrix[i - 1][j] + MMG.Gap, FSS.S_Matrix[i][j - 1] + MMG.Gap,
                                             FSS.S_Matrix[i - 1][j - 1] + MMG.Mismatch);
                }
            }
        }
    }
}

// Function for scoring and displaying the scoring matrix
void Scoring() {
    cout << "\n" << termcolor::on_cyan << "      Score Matrix :     " << termcolor::reset << "\n";
    for (int i = 0; i <= FSS.First_Seq.size(); i++) {
        for (int j = 0; j <= FSS.Second_Seq.size(); j++) {
            if (FSS.S_Matrix[i][j] >= 0) {
                cout << "  " << FSS.S_Matrix[i][j] << "      ";
            } else {
                cout << "  " << FSS.S_Matrix[i][j] << "      ";
            }
        }
        cout << "\n";
    }
    cout << termcolor::on_cyan << "\n The best score :" << termcolor::reset << " " << FSS.S_Matrix[FSS.First_Seq.size()][FSS.Second_Seq.size()] << "\n";
}

// Function for trace back and printing the aligned sequences
void traceBackANDPrint() {
    int i, j = 0;
    string val, val2;
    i = FSS.First_Seq.size();
    j = FSS.Second_Seq.size();

    while (i > 0 && j > 0) {
        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i - 1][j] + MMG.Gap) {
            PathT.push_back(FSS.S_Matrix[i - 1][j]);
            val = FSS.First_Seq[i - 1];
            val2 = FSS.Second_Seq[j - 1];
            Seq2.push_back(std::string(val2));
            Seq.push_back(std::string("-"));
            i = i - 1;
            j = j;
        }

        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i][j - 1] + MMG.Gap) {
            PathT.push_back(FSS.S_Matrix[i][j - 1]);
            val = FSS.First_Seq[i - 1];
            val2 = FSS.Second_Seq[j - 1];
            Seq.push_back(std::string(val));
            Seq2.push_back(std::string("-"));
            i = i;
            j = j - 1;
        }

        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i - 1][j - 1] + MMG.Match) {
            PathT.push_back(FSS.S_Matrix[i - 1][j - 1]);
            val = FSS.First_Seq[i - 1];
            val2 = FSS.Second_Seq[j - 1];
            Seq2.push_back(std::string(val2));
            Seq.push_back(std::string(val));
            i = i - 1;
            j = j - 1;
        }

        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i - 1][j - 1] + MMG.Mismatch) {
            PathT.push_back(FSS.S_Matrix[i - 1][j - 1]);
            val = FSS.First_Seq[i - 1];
            val2 = FSS.Second_Seq[j - 1];
            Seq2.push_back(std::string(val2));
            Seq.push_back(std::string(val));
            i = i - 1;
            j = j - 1;
        }
    }

    i = FSS.First_Seq.size();
    j = FSS.Second_Seq.size();

    int Counter = 0;

    while (i > 0 && j > 0) {
        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i - 1][j] + MMG.Gap) {
            if (FSS.S_Matrix[i][j] != PathT[Counter]) {
                PathT2.push_back(FSS.S_Matrix[i - 1][j]);
                val = FSS.First_Seq[i - 1];
                val2 = FSS.Second_Seq[j - 1];
                Seqq2.push_back(std::string(val2));
                Seqq.push_back(std::string("-"));
                i = i - 1;
                j = j;
            }
        }

        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i][j - 1] + MMG.Gap) {
            if (FSS.S_Matrix[i][j] != PathT[Counter]) {
                PathT2.push_back(FSS.S_Matrix[i][j - 1]);
                val = FSS.First_Seq[i - 1];
                val2 = FSS.Second_Seq[j - 1];
                Seqq.push_back(std::string(val));
                Seqq2.push_back(std::string("-"));
                i = i;
                j = j - 1;
            }
        }

        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i - 1][j - 1] + MMG.Match) {
            if (FSS.S_Matrix[i][j] != PathT[Counter]) {
                PathT2.push_back(FSS.S_Matrix[i - 1][j - 1]);
                val = FSS.First_Seq[i - 1];
                val2 = FSS.Second_Seq[j - 1];
                Seqq2.push_back(std::string(val2));
                Seqq.push_back(std::string(val));
                i = i - 1;
                j = j - 1;
            }
        }

        if (FSS.S_Matrix[i][j] == FSS.S_Matrix[i - 1][j - 1] + MMG.Mismatch) {
            if (FSS.S_Matrix[i][j] != PathT[Counter]) {
                PathT2.push_back(FSS.S_Matrix[i - 1][j - 1]);
                val = FSS.First_Seq[i - 1];
                val2 = FSS.Second_Seq[j - 1];
                Seqq2.push_back(std::string(val2));
                Seqq.push_back(std::string(val));
                i = i - 1;
                j = j - 1;
            }
        }

        Counter++;
    }

    // Print results
    cout << termcolor::on_cyan << "\n path 1 is :" << termcolor::reset;
    for (int i = 0; i < PathT.size(); i++) {
        cout << termcolor::green << " -> " << termcolor::reset << PathT[i];
    }
    cout << "\n";

    cout << termcolor::magenta << "\n Sequence 1 is :" << termcolor::reset;
    for (int i = 0; i < Seq.size(); i++) {
        cout << termcolor::green << " " << termcolor::reset << Seq[i];
    }
    cout << "\n";

    cout << termcolor::magenta << "\n Sequence 2 is :" << termcolor::reset;
    for (int i = 0; i < Seq2.size(); i++) {
        cout << termcolor::green << " " << termcolor::reset << Seq2[i];
    }
    cout << "\n";

    cout << termcolor::on_cyan << "\n path 2 is :" << termcolor::reset;
    for (int i = 0; i < PathT2.size(); i++) {
        cout << termcolor::green << " -> " << termcolor::reset << PathT2[i];
    }
    cout << "\n";

    cout << termcolor::magenta << "\n Sequence 1 is :" << termcolor::reset;
    for (int i = 0; i < Seqq.size(); i++) {
        cout << termcolor::green << " " << termcolor::reset << Seqq[i];
    }
    cout << "\n";

    cout << termcolor::magenta << "\n Sequence 2 is :" << termcolor::reset;
    for (int i = 0; i < Seqq2.size(); i++) {
        cout << termcolor::green << " " << termcolor::reset << Seqq2[i];
    }
}
