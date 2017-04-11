static
void FpLLL(mat_ZZ &B, int n)
{
    ofstream myfile;
    myfile.open ("mat.txt");
    myfile << B;
    myfile.close();

    char cmd[100];
    sprintf(cmd, "fplll -a lll mat.txt");

    FILE* pipe = popen(cmd, "r");

    string result = "";

    char buffer[100000];

    while (!feof(pipe))
    {
        if (fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }

    pclose(pipe);

    stringstream stream(result);
    stream >> B;
}


static
void FpBKZ(mat_ZZ &B, int n, int k)
{
    ofstream myfile;
    myfile.open ("mat.txt");
    myfile << B;
    myfile.close();

    char cmd[100];
    sprintf(cmd, "fplll -a bkz -b %d mat.txt",k);

    FILE* pipe = popen(cmd, "r");

    string result = "";

    char buffer[100000];

    while (!feof(pipe))
    {
        if (fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }

    pclose(pipe);

    stringstream stream(result);
    stream >> B;
}
