/* 3dcube.cpp
 *
 * Daniel Friedman
 * CS357
 * Fall 2012
 *
 * I affirm that I have adhered to the Honor Code in this assignment.
 */

#include <GL/glut.h>
#include <GL/glui.h>
#include <math.h>
#include <stdio.h>

#define true 1
#define false 0
#define WIN_SIZE 520
#define PI 3.1415926535

int main_window;

/* globals */
GLfloat SIZE = 1.0; /* size of the cube */
GLfloat av = 10.0, bv = 10.0, cv = 10.0; /* viewer position coordinates */
GLfloat lx = 0.0, ly = 0.0, lz = 0.0; /* viewer looking at coordinates */
/* clipping parameters */
GLfloat H = 5.0; /* distance between viewer and hither plane */
GLfloat D = 15.0; /* distance between viewer and picture plane */
GLfloat Y = 20.0; /* distance between viewer and yon plane */
GLfloat theta = 40.0; /* view angle */
/* viewport parameters */
GLfloat vr;
GLfloat vl;
GLfloat vt;
GLfloat vb;
/* window parameters */
GLfloat wr = 500.0;
GLfloat wl = 20.0;
GLfloat wt = 20.0;
GLfloat wb = 500.0;
/* matrices */
float V[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
float P[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
float W[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* math functions */


GLfloat toRadians(float deg)
{
    return deg*PI/180.0;
}

/* multiplies two 4x4 matricies (implemented as 1D arrays) */
void multiply(float a[16], float b[16], float result[16])
{
    int i, col, row; // the row and column of the product matrix
    for (row=0; row<4; row++) {
	for (col=0; col<4; col++) {
	    for (i=0; i<4; i++) { // get the dot product
                // using variables for debugging -- gdb
                // looks like result[(row*4)+col] is the bad value
                float x = a[(row*4)+i];
                float y = b[(i*4)+col];
                float prod = 0;
                prod = x * y;
		result[(row*4)+col] += prod;
	    }
	}
    }
}

/* multiplies a vector by a matrix, producing a vector */
void multiply2(float v[4], float m[16], float result[4])
{
    int row;
    int col;
    for (row = 0; row < 4; row++) {
        for (col = 0; col < 4; col++) {
            result[row] += v[col] * m[(col*4)+row];
        }
    }
}

/* computes V */
void view_transform()
{
    // compute az, bz, cz -- direction of view vector
    float az = lx - av;
    float bz = ly - bv;
    float cz = lz - cv;
    // compute ax, bx, cx -- take the cross product of <az,bz,cz> X up
    // can multiply scew matrix by up vector
    // always results in the following:
    float ax = bz;
    float bx = -az;
    float cx = 0;

    // compute r, R, and h
    float r = sqrt(ax*ax + bx*bx);
    float R = sqrt(ax*ax + bx*bx + cx*cx);
    float h = r * sqrt(az*az + bz*bz + cz*cz);

    float V1[16] = { 1, 0, 0, 0, 
                     0, 1, 0, 0, 
                     0, 0, 1, 0, 
                     -av, -bv, -cv, 1 };
    float V2[16] = { ax/r, (-bx)/r, 0, 0,
                     bx/r, ax/r, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1 };
    float V3[16] = { r/R, 0, (-cx)/R, 0,
                     0, 1, 0, 0,
                     cx/R, 0, r/R, 0,
                     0, 0, 0, 1 };
    float V4[16] = { 1, 0, 0, 0,
                     0, (cz*R)/h, (bz*ax-az*bx)/h, 0,
                     0, (az*bx-bz*ax)/h, (cz*R)/h, 0,
                     0, 0, 0, 1 };
    float V5[16] = { 1, 0, 0, 0,
                     0, 1, 0, 0, // 1 makes it look like Bob's example
                     0, 0, 1, 0,
                     0, 0, 0, 1 };
    float x1[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float x2[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    float x3[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    // reset V
    memset(V, 0, sizeof(int)*16);

    multiply(V1, V2, x1);
    multiply(x1, V3, x2);
    multiply(x2, V4, x3);
    multiply(x3, V5, V);
}

/* transforms viewer coordinates by perspective */
void persp_transform()
{
    float Ptemp[16] = { D, 0, 0, 0,
                        0, D, 0, 0,
                        0, 0, Y/(Y-H), 1,
                        0, 0, (-(H*Y))/(Y-H), 0
    };
    int i;
    for (i = 0; i < 16; i++) {
        P[i] = Ptemp[i];
    }
}

/* transforms viewer coordinates to window coordinates */
float* win_transform()
{
    vl = vb = -(D*tan(toRadians(theta)));
    vt = vr = D*tan(toRadians(theta));
    float Wtemp[16] = { (wr-wl)/(vr-vl), 0, 0, 0,
                        0, (wt-wb)/(vt-vb), 0, 0,
                        0, 0, 1, 0,
                        (wl*vr-vl*wr)/(vr-vl), (wb*vt-vb*wt)/(vt-vb), 0, 1
    };
    int i;
    for (i = 0; i < 16; i++) {
        W[i] = Wtemp[i];
    }
}

int clip(float a[4], float b[4]) // returns 0 if line discarded, 1 otherwise
{
    float *A;
    float *B;

    
    /* clip left */
    
    // if both a and b are to the left of wl, discard
    if (a[0] < wl && b[0] < wl) {
        return 1;
    }

    // figure out which point is to the left of which
    if (a[0] < b[0]) {
        A = a;
        B = b;
    } else {
        A = b;
        B = a;
    }
    
    // if a is to the left and b is to the right, find the intersection of the line with x=wl,
    if (A[0] < wl && B[0] >= wl) {
        float m = (B[1] - A[1]) / (B[0] - A[0]); // the slope
        float x = wl - A[0]; // horizontal distance from a to wl
        float y = m * x; // vertical distance from a to a'
        A[0] = wl; // change x coord to wl
        A[1] += y; // change y coord
    }

    
    /* clip right */
    // if both a and b are to the right of wr, discard
    if (a[0] > wr && b[0] > wr) {
        return 1;
    }

    // figure out which point is to the right of which
    if (a[0] > b[0]) {
        A = a;
        B = b;
    } else {
        A = b;
        B = a;
    }
    
    // if a is to the left and b is to the right, find the intersection of the line with x=wl,
    if (A[0] > wr && B[0] <= wr) {
        float m = (A[1] - B[1]) / (A[0] - B[0]); // the slope
        float x = A[0] - wr; // horizontal distance from a to wr
        float y = m * x; // vertical distance from a to a'
        A[0] = wr; // change x coord to wl
        A[1] -= y; // change y coord
    } 
    
    /* clip "top" -- openGL's coordinate system is ours upside-down */
    // if both a and b are above wr, discard
    if (a[1] < wt && b[1] < wt) {
        return 1;
    }

    // figure out which point is above which
    if (a[1] < b[1]) {
        A = a;
        B = b;
    } else {
        A = b;
        B = a;
    }
    
    // if a is above wt and b is below, find the intersection of the line with x=wt
    if (A[1] < wt && B[1] >= wt) {
        float m = (A[1] - B[1]) / (A[0] - B[0]); // the slope
        float y = wt - A[1]; // vertical distance from a to wt
        float x = y / m; // vertical distance from a to a'
        A[0] += x; // change x coord
        A[1] = wt; // change y coord to wt
    }

    /* clip "bottom" -- openGL's coordinate system is ours upside-down */
    // if both a and b are above wr, discard
    if (a[1] > wb && b[1] > wb) {
        return 1;
    }

    // figure out which point is above which
    if (a[1] > b[1]) {
        A = a;
        B = b;
    } else {
        A = b;
        B = a;
    }
    
    // if a is below wb and b is above, find the intersection of the line with x=wt
    if (A[1] > wb && B[1] <= wb) {
        float m = (A[1] - B[1]) / (A[0] - B[0]); // the slope
        float y = A[1] - wb; // vertical distance from a to wt
        float x = y / m; // vertical distance from a to a'
        A[0] += x; // change x coord
        A[1] = wb; // change y coord to wt
    }


    /* clip hither */
    // if both a and b are closer than H, discard
    if (a[2] < 0 && b[2] < 0) {
        return 1;
    }

    // figure out which point is closer
    if (a[2] < b[2]) {
        A = a;
        B = b;
    } else {
        A = b;
        B = a;
    }
    
    // if a is closer than H and b is farther from H, find the new point
    if (A[2] < 0 && B[2] >= 0) {
        float m1 = (A[1] - B[1]) / (A[2] - B[2]); // slope1
        float m2 = (A[2] - B[2]) / (A[0] - B[0]); // slope2
        float delta_z = 0 - A[2];
        float delta_y = m1 * delta_z;
        float delta_x = delta_z / m2;
        A[0] += delta_x; // change x coord
        A[1] += delta_y; // change y coord
        A[2] = 0;
    }


    /* clip yon */
    // if both a and b are further than Y, discard
    if (a[2] > 1 && b[2] > 1) {
        return 1;
    }

    // figure out which point is farther
    if (a[2] > b[2]) {
        A = a;
        B = b;
    } else {
        A = b;
        B = a;
    }
    
    // if a is further than Y and b is closer from Y, find the new point
    if (A[2] > 1 && B[2] <= 1) {
        float m1 = (A[1] - B[1]) / (A[2] - B[2]); // slope1
        float m2 = (A[0] - B[0]) / (A[2] - B[2]); // slope2
        float delta_z = A[2] - 1;
        float delta_y = m1 * delta_z;
        float delta_x = m2 * delta_z;
        A[0] -= delta_x; // change x coord
        A[1] -= delta_y; // change y coord
        A[2] = 1;
    }


    return 0;
}



/* drawing functions */

/* draws a line from specified workd coordinates */
void drawline(float *p1, float *p2)
{
    // compute the screen coordinate for p1
    float vV[4] = {0, 0, 0, 0};
    multiply2(p1, V, vV);
    float vVP[4] = {0, 0, 0, 0};
    multiply2(vV, P, vVP);
    float prod1[4] = {0, 0, 0, 0};
    multiply2(vVP, W, prod1);
    // homogenize result of <p1>*VPW
    int i;
    for (i = 0; i < 4; i++) {
        prod1[i] /= prod1[3];
    }

    // compute the screen coordinate for p2
    float vV2[4] = {0, 0, 0, 0};
    multiply2(p2, V, vV2);
    float vVP2[4] = {0, 0, 0, 0};
    multiply2(vV2, P, vVP2);
    float prod2[4] = {0, 0, 0, 0};
    multiply2(vVP2, W, prod2);
    // homogenize result of <p2>*VPW
    for (i = 0; i < 4; i++) {
        prod2[i] /= prod2[3];
    }

    // clip
    if(!clip(prod1, prod2)) {
        glBegin(GL_LINES);
            glVertex2f(prod1[0], prod1[1]);
            glVertex2f(prod2[0], prod2[1]);
        glEnd();
    }
}




/* program functions */

/* sets up the window */
void setupViewport(int w, int h)
{
    glViewport(0, 0, w, h); 
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluOrtho2D(0.0, 520*w/WIN_SIZE, 0.0, 520*h/WIN_SIZE );
}

void reshape(int w, int h)
{
    setupViewport(w, h);
    glutPostWindowRedisplay(main_window);
}


/* draws a viewport onto the window so we can see clipping */
void drawPort()
{
    glColor3f(0.0, 0.4, 0.4); // light magenta

    glBegin(GL_POLYGON);
        glVertex2f(20, 20);
        glVertex2f(WIN_SIZE - 20, 20);
        glVertex2f(WIN_SIZE - 20, WIN_SIZE - 20);
        glVertex2f(20, WIN_SIZE - 20);
    glEnd();
    glutSwapBuffers();
}

void init()
{
    glClearColor(1.0, 1.0, 1.0, 0.0); // set the background color to white
    glColor3f(1.0, 0.0, 0.0);

    // initialize V, P, and W
    view_transform();
    persp_transform();
    win_transform();

    // prints out the values of V, P, and W for debugging
    //printf("V matrix:\n");
    //for (int i = 0; i < 16; i++) {
    //    if (i%4 == 0) {
    //        printf("\n");
    //    }
    //    printf("[%f],\t", V[i]);
    //}
    //printf("\n");
    //printf("P matrix:\n");
    //for (int i = 0; i < 16; i++) {
    //    if (i%4 == 0) {
    //        printf("\n");
    //    }
    //    printf("[%f],\t", P[i]);
    //}
    //printf("\n");
    //printf("W matrix:\n");
    //for (int i = 0; i < 16; i++) {
    //    if (i%4 == 0) {
    //        printf("\n");
    //    }
    //    printf("[%f],\t", W[i]);
    //}
    //printf("\n");

    setupViewport(520, 520); // setup our initial window
}

void display()
{
    glutSetWindow(main_window); // chooses the correct window
    glClear(GL_COLOR_BUFFER_BIT); // clears the buffer

    drawPort();

    // compute points
    float p1[4] = {0, 0, 0, 1};
    float p2[4] = {0, 0, 1*SIZE, 1};
    float p3[4] = {0, 1*SIZE, 0, 1};
    float p4[4] = {0, 1*SIZE, 1*SIZE, 1};
    float p5[4] = {1*SIZE, 0, 0, 1};
    float p6[4] = {1*SIZE, 0, 1*SIZE, 1};
    float p7[4] = {1*SIZE, 1*SIZE, 0, 1};
    float p8[4] = {1*SIZE, 1*SIZE, 1*SIZE, 1};


    // cube edges
    glColor3f(0.0, 0.0, 1.0); // blue
    drawline(p1, p2);
    drawline(p1, p3);
    drawline(p2, p4);
    drawline(p3, p4);
    glColor3f(1.0, 0.0, 0.0); // red
    drawline(p5, p6);
    drawline(p6, p8);
    drawline(p5, p7);
    drawline(p7, p8);
    glColor3f(0.0, 1.0, 0.0); // green
    drawline(p1, p5);
    drawline(p3, p7);
    drawline(p4, p8);
    drawline(p2, p6);
    drawline(p5, p8); // diagonal

    glutSwapBuffers(); // draws the new buffer
}




/* callback functions */

/* whenever viewer controls are changed */
void viewer_callback(int ID)
{
    view_transform();
    display(); // draw the cube
}

/* when theta is changed */
void theta_callback(int ID)
{
    win_transform();
    display(); // draw the cube
}

/* when H, D, or Y are changed */
void clipping_callback(int ID)
{
    persp_transform();
    display();
}

/* callback for when window is resized */
void window_callback(int ID)
{
    win_transform();
    display(); // draw the cube
}


/* test functions */
void test1()
{
    av = 10.0;
    bv = 10.0;
    cv = 10.0;
    lx = ly = lz = 0.0;

    view_transform();
    persp_transform();
    win_transform();

    float point[4] = {2.0, 5.0, 1.0, 1.0};

    float vV[4] = {0, 0, 0, 0};
    multiply2(point, V, vV);
    printf("(2, 5, 1)V = [%f, %f, %f, %f]\n", vV[0], vV[1], vV[2], vV[3]);
    float vVP[4] = {0, 0, 0, 0};
    multiply2(vV, P, vVP);
    printf("(2, 5, 1)VP = [%f, %f, %f, %f]\n", vVP[0], vVP[1], vVP[2], vVP[3]);
    float prod1[4] = {0, 0, 0, 0};
    multiply2(vVP, W, prod1);

    printf("(2, 5, 1)VPW = [%f, %f, %f, %f]\n", prod1[0], prod1[1], prod1[2], prod1[3]);
    // homogenize result of <p1>*VPW
    int i;
    for (i = 0; i < 4; i++) {
        prod1[i] /= prod1[3];
    }

    printf("(2, 5, 1)VPW<homogenize> = [%f, %f, %f, %f]\n", prod1[0], prod1[1], prod1[2], prod1[3]);
}




int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // use double-buffering, RGB
    glutInitWindowSize(WIN_SIZE, WIN_SIZE);
    glutInitWindowPosition(50, 50);
    main_window = glutCreateWindow("3D Cube");

    init(); // call our init function
    glutDisplayFunc(display); // register display function
    glutReshapeFunc(reshape);

    // create the control panel, initialize some buttons
    GLUI *control_panel = GLUI_Master.create_glui("Control box");
    new GLUI_StaticText(control_panel, "3D Cube");
    new GLUI_Separator(control_panel);
    new GLUI_Button(control_panel, "Quit", 0, (GLUI_Update_CB) exit);
    new GLUI_Column(control_panel, true);

    // Size panel
    GLUI_Spinner *size = new GLUI_Spinner(control_panel, "SIZE", GLUI_SPINNER_FLOAT, &SIZE, 0, viewer_callback);
    size->set_float_limits(1.0, 20.0);
    new GLUI_Column(control_panel, true);
    // eye position / looking at panel
    // eye position rollout
    GLUI_Rollout *eye_pos_rollout = new GLUI_Rollout(control_panel, "Eye Position", false);
    GLUI_Spinner *x_pos= new GLUI_Spinner(eye_pos_rollout, "X", GLUI_SPINNER_FLOAT, &av, 0, viewer_callback);
    GLUI_Spinner *y_pos= new GLUI_Spinner(eye_pos_rollout, "Y", GLUI_SPINNER_FLOAT, &bv, 0, viewer_callback);
    GLUI_Spinner *z_pos= new GLUI_Spinner(eye_pos_rollout, "Z", GLUI_SPINNER_FLOAT, &cv, 0, viewer_callback);
    // looking at rollout
    GLUI_Rollout *look_at_rollout = new GLUI_Rollout(control_panel, "Looking At", false);
    GLUI_Spinner *x_look= new GLUI_Spinner(look_at_rollout, "X", GLUI_SPINNER_FLOAT, &lx, 0, viewer_callback);
    GLUI_Spinner *y_look= new GLUI_Spinner(look_at_rollout, "Y", GLUI_SPINNER_FLOAT, &ly, 0, viewer_callback);
    GLUI_Spinner *z_look= new GLUI_Spinner(look_at_rollout, "Z", GLUI_SPINNER_FLOAT, &lz, 0, viewer_callback);

    new GLUI_Column(control_panel, true);
    // clipping param panel
    GLUI_Rollout *clipping_rollout = new GLUI_Rollout(control_panel, "Clipping Parameters plane", false);
    GLUI_Spinner *dist_hither = new GLUI_Spinner(clipping_rollout, "Distance to Hither plane", GLUI_SPINNER_FLOAT, &H, 0, clipping_callback);
    dist_hither->set_float_limits(1.0, 20.0, GLUI_LIMIT_CLAMP);
    GLUI_Spinner *dist_picture = new GLUI_Spinner(clipping_rollout, "Distance to Picture plane", GLUI_SPINNER_FLOAT, &D, 0, clipping_callback);
    dist_picture->set_float_limits(5.0, 25.0, GLUI_LIMIT_CLAMP);
    GLUI_Spinner *dist_yon = new GLUI_Spinner(clipping_rollout, "Distance to Yon plane", GLUI_SPINNER_FLOAT, &Y, 0, clipping_callback);
    dist_yon->set_float_limits(10.0, 30.0, GLUI_LIMIT_CLAMP);


    // theta spinner
    GLUI_Spinner *theta_spinner= new GLUI_Spinner(control_panel, "THETA", GLUI_SPINNER_FLOAT, &theta, 0, theta_callback);
    theta_spinner->set_float_limits(30.0, 75.0, GLUI_LIMIT_CLAMP);

    control_panel->set_main_gfx_window(main_window);

    /* call the test function (for... testing) */
    //test1();
    
    //GLUI_Master.set_glutIdleFunc(spin);
    glutMainLoop();

    return EXIT_SUCCESS;
}
