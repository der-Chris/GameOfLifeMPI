#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <mpi.h>
#include <time.h>

#define calcIndex(width, x,y)  ((y + 1)*(width + 2) + (x + 1))

void show(unsigned* currentfield, int w, int h) {
    printf("\033[H");
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
        printf("\033[E");
    }
    fflush(stdout);
}


float convert2BigEndian( const float inFloat )
{
    float retVal;
    char *floatToConvert = ( char* ) & inFloat;
    char *returnFloat    = ( char* ) & retVal;

    // swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];

    return retVal;
}

void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix, MPI_Comm comm, int rank, int size, int *myCoords) {
    char name[1024] = "\0";
    sprintf(name, "out/%s_%d_%d.vtk", prefix, rank, t);
    FILE* outfile = fopen(name, "w");

    int tempy = (h - 1) * myCoords[0];
    int tempx = (w - 1) * myCoords[1];

    /*Write vtk header */
    fprintf(outfile,"# vtk DataFile Version 3.0\n");
    fprintf(outfile,"frame %d\n", t);
    fprintf(outfile,"BINARY\n");
    fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
    fprintf(outfile,"DIMENSIONS %d %d %d \n", w, h, 1);
    fprintf(outfile,"SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO
    fprintf(outfile,"ORIGIN %d %d 0\n", tempx, tempy);
    fprintf(outfile,"POINT_DATA %d\n", h*w);
    fprintf(outfile,"SCALARS data float 1\n");
    fprintf(outfile,"LOOKUP_TABLE default\n");

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float value = currentfield[calcIndex(w, x,y)]; // != 0.0 ? 1.0:0.0;
            value = convert2BigEndian(value);
            fwrite(&value, 1, sizeof(float), outfile);
        }
    }
    fclose(outfile);
}

void calculateNeighbourProcesses(int *myCoords, MPI_Comm comm, int *rightProcess, int *leftProcess)
{
    int *queryCoords = calloc(2, sizeof(unsigned));

    queryCoords[0] = myCoords[0];
    queryCoords[1] = myCoords[1] + 1;
    MPI_Cart_rank(comm, queryCoords, rightProcess);

    queryCoords[0] = myCoords[0];
    queryCoords[1] = myCoords[1] - 1;
    MPI_Cart_rank(comm, queryCoords, leftProcess);

    free(queryCoords);
}

void calculateOwnBorders(unsigned *currentfield, int myW, int myH, unsigned *sendRightBorder, unsigned *sendLeftBorder) {
    for (int i = 0; i < myH; i = i + 1)
    {
        sendRightBorder[i] = currentfield[calcIndex(myW, myW - 1, i)];
        sendLeftBorder[i] = currentfield[calcIndex(myW, 0, i)];
    }
}

void calculateNeighbourBorders(unsigned *currentfield, int w, int h, MPI_Comm comm, unsigned *sendRightBorder, unsigned *sendLeftBorder, unsigned *recvRightBorder, unsigned *recvLeftBorder, int *rightProcess, int *leftProcess) {

    MPI_Status status;

    //Send Upper Border -> receive lower Border
    MPI_Datatype sendtype = MPI_UNSIGNED;
    int sendtag = 1;
    int recvtag = 1;
    MPI_Datatype recvtype = MPI_UNSIGNED;
    //Send Lower Border -> receive Upper Border

    //Send Right Border -> receive Left Border
    int sendcount = h;
    int dest = rightProcess[0];
    int recvcount = h;
    int source = leftProcess[0];
    MPI_Sendrecv(sendRightBorder, sendcount, sendtype, dest, sendtag, recvLeftBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Left Border -> receive Right Border
    dest = leftProcess[0];
    source = rightProcess[0];
    MPI_Sendrecv(sendLeftBorder, sendcount, sendtype, dest, sendtag, recvRightBorder, recvcount, recvtype, source, recvtag, comm, &status);
}

int calcAlive(unsigned *currentfield, unsigned *newfield, int myW, int myH)
{
    int changes = 0;
    for (int y = 0; y < myH; y++) {
        for (int x = 0; x < myW; x++) {
            int alive = 0;
            //for schleifen ersetzen durch 8 statische berechnungen?
            for(int y1 = y-1; y1 <= y+1; y1++) {
                for(int x1 = x-1; x1 <= x+1; x1++) {
                    if(currentfield[calcIndex(myW, x1, y1)]) {
                        alive++;
                    }
                }
            }

            //remove inner cell value
            alive -= currentfield[calcIndex(myW, x,y)];
            newfield[calcIndex(myW, x,y)] = (alive == 3 || (alive == 2 && currentfield[calcIndex(myW, x,y)]));
            if(currentfield[calcIndex(myW, x,y)] != newfield[calcIndex(myW, x,y)])
            {
                changes = 1;
            }
        }
    }
    return changes;
}

void mergeBorders(unsigned *currentfield, int myW, int myH, unsigned *recvRightBorder, unsigned *recvLeftBorder) {
    for (int i = 0; i < myH; i = i + 1)
    {
        currentfield[calcIndex(myW, myW, i)] = recvRightBorder[i];
        currentfield[calcIndex(myW, -1, i)] = recvLeftBorder[i];
    }
}

int evolve(unsigned *currentfield, unsigned *newfield, MPI_Status status, MPI_Comm comm, int rank, int size, int *myCoords, int w, int h) {
    int changes = 0;

    unsigned *sendRightBorder = calloc(h, sizeof(unsigned));
    unsigned *sendLeftBorder = calloc(h, sizeof(unsigned));

    unsigned *recvRightBorder = calloc(h, sizeof(unsigned));
    unsigned *recvLeftBorder = calloc(h, sizeof(unsigned));

    int *rightProcess = calloc (1, sizeof(int));
    int *leftProcess = calloc (1, sizeof(int));

    calculateNeighbourProcesses(myCoords, comm, rightProcess, leftProcess);

    calculateOwnBorders(currentfield, w, h, sendRightBorder, sendLeftBorder);

    calculateNeighbourBorders(currentfield, w, h, comm, sendRightBorder, sendLeftBorder, recvRightBorder, recvLeftBorder, rightProcess, leftProcess);

    mergeBorders(currentfield, w, h, recvRightBorder, recvLeftBorder);

    changes = calcAlive(currentfield, newfield, w, h);

    free(sendRightBorder);
    free(sendLeftBorder);

    free(recvRightBorder);
    free(recvLeftBorder);

    free(rightProcess);
    free(leftProcess);

    //TODO if changes == 0, the time loop will not run!
    return changes;
}

void filling(unsigned* currentfield, int w, int h, int rank) {
    clock_t t = time(NULL);
    srand((int) t * rank);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            currentfield[calcIndex(w, x, y)] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
        }
    }
}

void game(int w, int h, int timesteps, int rank, int size, MPI_Comm comm, MPI_Status status) {
    w = w / size;
    unsigned *currentfield = calloc((w + 2) * (h + 2), sizeof(unsigned));
    unsigned *newfield = calloc((w + 2) * (h + 2), sizeof(unsigned));

    MPI_Comm newComm;
    int ndims = 2;
    int *dims = calloc(ndims, sizeof(unsigned));
    dims[0] = 1;
    dims[1] = size;
    int *periods = calloc(ndims, sizeof(int));
    periods[0] = 1;
    periods[1] = 1;
    int reorder = 0;

    MPI_Cart_create(comm, ndims, dims, periods, reorder, &newComm);
    int *myCoords = calloc(ndims, sizeof(int));
    MPI_Cart_coords(newComm, rank, ndims, myCoords);

    //printf("DEBUG: Hello, my ID is [%d],  my right neighbours ID is [%d] my left neighbours ID is [%d]\n", rank, rightProcess[0], leftProcess[0]);
    printf("DEBUG: Hello, my ID is [%d],  my field is %dX%d\n", rank, w, h);

    filling(currentfield, w, h, rank);

    for (int t = 0; t < timesteps; t++) {
        int t = 0;
        writeVTK(currentfield, w, h, t, "output", newComm, rank, size, myCoords);
        int changes = evolve(currentfield, newfield, status, newComm, rank, size, myCoords, w, h);
        int allchanges = changes;
        // printf("Changes of Thread %d is %d\n", rank, changes);
        MPI_Allreduce(&changes, &allchanges, 1, MPI_INT, MPI_SUM, comm);
        // printf("Changes of all Threads is %d\n", allchanges);
        if (allchanges == 0) {
            sleep(3);
            break;
        }

        //SWAP
        unsigned *temp = currentfield;
        currentfield = newfield;
        newfield = temp;
    }

    free(currentfield);
    free(newfield);
}

int main(int c, char **v) {
    int rank, size;
    MPI_Status status;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&c, &v);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);


    int w = 0, h = 0, timesteps = 10;
    if (c > 1) w = atoi(v[1]); ///< read width
    if (c > 2) h = atoi(v[2]); ///< read height
    if (c > 3) timesteps = atoi(v[3]);
    if (w <= 0) w = 32; ///< default width
    if (h <= 0) h = 32; ///< default height
    game(w, h, timesteps, rank, size, comm, status);
    MPI_Finalize();
    return 0;
}
