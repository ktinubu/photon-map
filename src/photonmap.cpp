// Source file for the scene viewer program

// https://stackoverflow.com/questions/12702561/iterate-through-a-c-vector-using-a-for-loop


// Include files 

#include "R3Graphics/R3Graphics.h"
#include "fglut/fglut.h"
#include "render.h"
#include <iostream>
#include <vector>
#include <initializer_list>
#include <string>





// Program variables

static char *input_scene_name = NULL;
static char *output_image_name = NULL;
static char *screenshot_image_name = NULL;
static int print_verbose = 0;
static int max_bounces = -1;
static RNScalar termination_rate = 0.05; // rate at which photons get terminated
static RNScalar camera_index_of_refraction = 1.0;

// GLUT variables 
static int GLUTwindow = 1;
static int GLUTwindow_height = 800;
static int GLUTwindow_width = 800;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmouse_drag = 0;
static int GLUTmodifiers = 0;

// Application variables
static R3Viewer *viewer = NULL;
static R3Scene *scene = NULL;
static R3Point center(0, 0, 0);

// Photon mapping variables
static int render_image_width = 200;
static int render_image_height = 200;
static int num_photons = 500000;
static int num_caustics = 1000000;
static int num_samples = 20;
static RNScalar tone_map_const = 0.3;
static RNScalar general_search_range = 0.07; // as a proprtion of radius of bounding box of scene
static RNScalar caustic_search_range = 0.1; // as a proprtion of radius of bounding box of scene
static int num_photon_estimate = 150;


static RNArray<Photon *> photon_list;
static R3Kdtree<Photon *> *photon_map;
static R3Kdtree<Photon *> *caustic_map;
static RNArray<Photon *> caustic_list;

// Display variables

static int show_shapes = 1;
static int show_camera = 1;
static int show_lights = 1;
static int show_bboxes = 0;
static int show_rays = 0;
static int show_photons = 0;
static int show_caustic_photons = 0;
static int show_frame_rate = 0;
static int show_near_photons = 0;

////////////////////////////////////////////////////////////////////////
// Draw functions
////////////////////////////////////////////////////////////////////////

static void 
LoadLights(R3Scene *scene)
{
  // Load ambient light
  static GLfloat ambient[4];
  ambient[0] = scene->Ambient().R();
  ambient[1] = scene->Ambient().G();
  ambient[2] = scene->Ambient().B();
  ambient[3] = 1;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  // Load scene lights
  for (int i = 0; i < scene->NLights(); i++) {
    R3Light *light = scene->Light(i);
    light->Draw(i);
  }
}



#if 0

static void 
DrawText(const R3Point& p, const char *s)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}
  
#endif



static void 
DrawText(const R2Point& p, const char *s)
{
  // Draw text string s and position p
  R3Ray ray = viewer->WorldRay((int) p[0], (int) p[1]);
  R3Point position = ray.Point(2 * viewer->Camera().Near());
  glRasterPos3d(position[0], position[1], position[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}



static void 
DrawCamera(R3Scene *scene)
{
  // Draw view frustum
  const R3Camera& camera = scene->Camera();
  R3Point eye = camera.Origin();
  R3Vector towards = camera.Towards();
  R3Vector up = camera.Up();
  R3Vector right = camera.Right();
  RNAngle xfov = camera.XFOV();
  RNAngle yfov = camera.YFOV();
  double radius = scene->BBox().DiagonalRadius();
  R3Point org = eye + towards * radius;
  R3Vector dx = right * radius * tan(xfov);
  R3Vector dy = up * radius * tan(yfov);
  R3Point ur = org + dx + dy;
  R3Point lr = org + dx - dy;
  R3Point ul = org - dx + dy;
  R3Point ll = org - dx - dy;
  glBegin(GL_LINE_LOOP);
  glVertex3d(ur[0], ur[1], ur[2]);
  glVertex3d(ul[0], ul[1], ul[2]);
  glVertex3d(ll[0], ll[1], ll[2]);
  glVertex3d(lr[0], lr[1], lr[2]);
  glVertex3d(ur[0], ur[1], ur[2]);
  glVertex3d(eye[0], eye[1], eye[2]);
  glVertex3d(lr[0], lr[1], lr[2]);
  glVertex3d(ll[0], ll[1], ll[2]);
  glVertex3d(eye[0], eye[1], eye[2]);
  glVertex3d(ul[0], ul[1], ul[2]);
  glEnd();
}


static void 
DrawPhotons(R3Scene *scene, RNArray<Photon *> photon_list, R3Kdtree<Photon *> *photon_map)
{
  // Draw all lights
  double radius = scene->BBox().DiagonalRadius();
  R3Point camera_pos = scene->Camera().Origin();
  for (int i = 0; i < photon_list.NEntries(); i++) {
    Photon *photon = photon_list[i];
    int ab = 100000;
    glColor3d(photon->power[0] * ab, photon->power[1] * ab, photon->power[2] * ab);
    R3Sphere(photon->position, 0.005 * radius).Draw();
    R3Span(photon->position, photon->source).Draw();
    R3Span(photon->position, photon->position + 0.05 * radius * photon->normal).Draw();
    
    // R3Vector mirrored = photon->direction - (2 * photon->direction.Dot(photon->normal) * photon->normal);
   
    glColor3d(0.0, 1.0, 0.0);
  }
}

static void 
DrawNearestPhotons(R3Scene *scene, RNArray<Photon *> photon_list, R3Kdtree<Photon *> *photon_map)
{
  // Draw all lights
  double radius = scene->BBox().DiagonalRadius();
  R3Point camera_pos = scene->Camera().Origin();

  for (int i = 0; i < 5; i++) {
    Photon *source_photon = photon_list[i];
    R3Span(camera_pos, source_photon->position).Draw();
    RNArray<Photon *> nearby_photons;
    photon_map->FindClosest(source_photon, RNScalar(0), radius, 500, nearby_photons);
    if (nearby_photons.NEntries() == 0) {
      continue;
    }

    for (int p = 0; p < nearby_photons.NEntries(); p++) {
      Photon *near_photon = nearby_photons[p];
      int ab = 100000;
      glColor3d(near_photon->power[0] * ab, near_photon->power[1] * ab, near_photon->power[2] * ab);
      R3Sphere(near_photon->position, 0.005 * radius).Draw();
    }
    glColor3d(0.0, 1.0, 0.0);
  }
}
static void 
DrawLights(R3Scene *scene)
{
  // Draw all lights
  double radius = scene->BBox().DiagonalRadius();
  for (int i = 0; i < scene->NLights(); i++) {
    R3Light *light = scene->Light(i);
    RNLoadRgb(light->Color());
    if (light->ClassID() == R3DirectionalLight::CLASS_ID()) {
      R3DirectionalLight *directional_light = (R3DirectionalLight *) light;
      R3Vector direction = directional_light->Direction();

      // Draw direction vector
      glBegin(GL_LINES);
      R3Point centroid = scene->BBox().Centroid();
      R3LoadPoint(centroid - radius * direction);
      R3LoadPoint(centroid - 1.25 * radius * direction);
      glEnd();
    }
    else if (light->ClassID() == R3PointLight::CLASS_ID()) {
      // Draw sphere at point light position
      R3PointLight *point_light = (R3PointLight *) light;
      R3Point position = point_light->Position();

     // Draw sphere at light position 
       R3Sphere(position, 0.1 * radius).Draw();
    }
    else if (light->ClassID() == R3SpotLight::CLASS_ID()) {
      R3SpotLight *spot_light = (R3SpotLight *) light;
      R3Point position = spot_light->Position();
      R3Vector direction = spot_light->Direction();

      // Draw sphere at light position 
      R3Sphere(position, 0.1 * radius).Draw();
  
      // Draw direction vector
      glBegin(GL_LINES);
      R3LoadPoint(position);
      R3LoadPoint(position + 0.25 * radius * direction);
      glEnd();
    }
    else {
      fprintf(stderr, "Unrecognized light type: %d\n", light->ClassID());
      return;
    }
  }
}



static void 
DrawShapes(R3Scene *scene, R3SceneNode *node, RNFlags draw_flags = R3_DEFAULT_DRAW_FLAGS)
{
  // Push transformation
  node->Transformation().Push();

  // Draw elements 
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);
    element->Draw(draw_flags);
  }

  // Draw children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    DrawShapes(scene, child, draw_flags);
  }

  // Pop transformation
  node->Transformation().Pop();
}



static void 
DrawBBoxes(R3Scene *scene, R3SceneNode *node)
{
  // Draw node bounding box
  node->BBox().Outline();

  // Push transformation
  node->Transformation().Push();

  // Draw children bboxes
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    DrawBBoxes(scene, child);
  }

  // Pop transformation
  node->Transformation().Pop();
}



static void 
DrawRays(R3Scene *scene)
{
  // Ray intersection variables
  R3SceneNode *node;
  R3SceneElement *element;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // Ray generation variables
  int istep = scene->Viewport().Width() / 20;
  int jstep = scene->Viewport().Height() / 20;
  if (istep == 0) istep = 1;
  if (jstep == 0) jstep = 1;

  // Ray drawing variables
  double radius = 0.025 * scene->BBox().DiagonalRadius();

  // Draw intersection point and normal for some rays
  for (int i = istep/2; i < scene->Viewport().Width(); i += istep) {
    for (int j = jstep/2; j < scene->Viewport().Height(); j += jstep) {
      R3Ray ray = scene->Viewer().WorldRay(i, j);
      if (scene->Intersects(ray, &node, &element, &shape, &point, &normal, &t)) {
        R3Sphere(point, radius).Draw();
        R3Span(point, point + 2 * radius * normal).Draw();
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Glut user interface functions
////////////////////////////////////////////////////////////////////////

void GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



void GLUTRedraw(void)
{
  // Check scene
  if (!scene) return;

  // Set viewing transformation
  viewer->Camera().Load();

  // Clear window 
  RNRgb background = scene->Background();
  glClearColor(background.R(), background.G(), background.B(), 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Load lights
  LoadLights(scene);

  // Draw camera
  if (show_camera) {
    glDisable(GL_LIGHTING);
    glColor3d(1.0, 1.0, 1.0);
    glLineWidth(5);
    DrawCamera(scene);
    glLineWidth(1);
  }

  // Draw lights
  if (show_lights) {
    glDisable(GL_LIGHTING);
    glColor3d(1.0, 1.0, 1.0);
    glLineWidth(5);
    DrawLights(scene);
    glLineWidth(1);
  }

  // Draw rays
  if (show_rays) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    glLineWidth(3);
    DrawRays(scene);
    glLineWidth(1);
  }

  // Draw photons
  if (show_photons) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    glLineWidth(3);
    DrawPhotons(scene, photon_list, photon_map);
    glLineWidth(1);
  }

  // Draw photons
  if (show_caustic_photons) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    glLineWidth(3);
    DrawPhotons(scene, caustic_list, caustic_map);
    glLineWidth(1);
  }

  // Draw photons
  if (show_near_photons) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    glLineWidth(3);
    DrawNearestPhotons(scene, photon_list, photon_map);
    DrawNearestPhotons(scene, caustic_list, caustic_map);
    glLineWidth(1);
  }

  // Draw scene nodes
  if (show_shapes) {
    glEnable(GL_LIGHTING);
    R3null_material.Draw();
    DrawShapes(scene, scene->Root());
    R3null_material.Draw();
  }

  // Draw bboxes
  if (show_bboxes) {
    glDisable(GL_LIGHTING);
    glColor3d(1.0, 0.0, 0.0);
    DrawBBoxes(scene, scene->Root());
  }

  // Draw frame time
  if (show_frame_rate) {
    char buffer[128];
    static RNTime last_time;
    double frame_time = last_time.Elapsed();
    last_time.Read();
    if ((frame_time > 0) && (frame_time < 10)) {
      glDisable(GL_LIGHTING);
      glColor3d(1.0, 1.0, 1.0);
      sprintf(buffer, "%.1f fps", 1.0 / frame_time);
      DrawText(R2Point(100, 100), buffer);
    }
  }  

  // Capture screenshot image 
  if (screenshot_image_name) {
    if (print_verbose) printf("Creating image %s\n", screenshot_image_name);
    R2Image image(GLUTwindow_width, GLUTwindow_height, 3);
    image.Capture();
    image.Write(screenshot_image_name);
    screenshot_image_name = NULL;
  }

  // Swap buffers 
  glutSwapBuffers();
}    



void GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  viewer->ResizeViewport(0, 0, w, h);

  // Resize scene viewport
  scene->SetViewport(viewer->Viewport());

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}



void GLUTMotion(int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute mouse movement
  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];
  
  // Update mouse drag
  GLUTmouse_drag += dx*dx + dy*dy;

  // World in hand navigation 
  if (GLUTbutton[0]) viewer->RotateWorld(1.0, center, x, y, dx, dy);
  else if (GLUTbutton[1]) viewer->ScaleWorld(1.0, center, x, y, dx, dy);
  else if (GLUTbutton[2]) viewer->TranslateWorld(1.0, center, x, y, dx, dy);
  if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



void GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Mouse is going down
  if (state == GLUT_DOWN) {
    // Reset mouse drag
    GLUTmouse_drag = 0;
  }
  else {
    // Check for double click  
    static RNBoolean double_click = FALSE;
    static RNTime last_mouse_up_time;
    double_click = (!double_click) && (last_mouse_up_time.Elapsed() < 0.4);
    last_mouse_up_time.Read();

    // Check for click (rather than drag)
    if (GLUTmouse_drag < 100) {
      // Check for double click
      if (double_click) {
        // Set viewing center point 
        R3Ray ray = viewer->WorldRay(x, y);
        R3Point intersection_point;
        if (scene->Intersects(ray, NULL, NULL, NULL, &intersection_point)) {
          center = intersection_point;
        }
      }
    }
  }

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Redraw
  glutPostRedisplay();
}



void GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}



void GLUTKeyboard(unsigned char key, int x, int y)
{
  // Process keyboard button event 
  switch (key) {
  case '~': {
    // Dump screen shot to file iX.jpg
    static char buffer[64];
    static int image_count = 1;
    sprintf(buffer, "i%d.jpg", image_count++);
    screenshot_image_name = buffer;
    break; }

  case 'B':
  case 'b':
    show_bboxes = !show_bboxes;
    break;

  case 'C':
  case 'c':
    show_camera = !show_camera;
    break;

  case 'L':
  case 'l':
    show_lights = !show_lights;
    break;

  case 'R':
  case 'r':
    show_rays = !show_rays;
    break;

  case 'S':
  case 's':
    show_shapes = !show_shapes;
    break;

  case 'P':
  case 'p':
    show_photons = !show_photons;
    break;

  case 'Q':
  case 'q':
    show_caustic_photons = !show_caustic_photons;
    break;

  case 'T':
  case 't':
    show_frame_rate = !show_frame_rate;
    break;

  case ' ':
    viewer->SetCamera(scene->Camera());
    break;

  case 27: // ESCAPE
    GLUTStop();
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = GLUTwindow_height - y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();  
}




void GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
  GLUTwindow = glutCreateWindow("Property Viewer");

  // Initialize lighting
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING); 

  // Initialize headlight
  // static GLfloat light0_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  // static GLfloat light0_position[] = { 0.0, 0.0, 1.0, 0.0 };
  // glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  // glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  // glEnable(GL_LIGHT0);

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);
}



void GLUTMainLoop(void)
{
  // Initialize viewing center
  if (scene) center = scene->BBox().Centroid();

  // Run main loop -- never returns 
  glutMainLoop();
}


 
////////////////////////////////////////////////////////////////////////
// Input/output
////////////////////////////////////////////////////////////////////////

static R3Scene *
ReadScene(char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene for %s\n", filename);
    return NULL;
  }

  // Read scene from file
  if (!scene->ReadFile(filename)) {
    delete scene;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Nodes = %d\n", scene->NNodes());
    printf("  # Lights = %d\n", scene->NLights());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



static int
WriteImage(R2Image *image, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write image to file
  if (!image->Write(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote image to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Width = %d\n", image->Width());
    printf("  Height = %d\n", image->Height());
    fflush(stdout);
  }

  // Return success
  return 1;
}  



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse argumentstone_map_const
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) {
        print_verbose = 1; 
      } else if (!strcmp(*argv, "-resolution")) { 
        argc--; argv++; render_image_width = atoi(*argv); 
        argc--; argv++; render_image_height = atoi(*argv); 
      } else if (!strcmp(*argv, "-num_samples")) { 
        argc--; argv++; num_samples = atoi(*argv); 
      } else if (!strcmp(*argv, "-general_search_range")) { 
        argc--; argv++; general_search_range = atof(*argv); 
      } else if (!strcmp(*argv, "-caustic_search_range")) { 
        argc--; argv++; caustic_search_range = atof(*argv); 
      } else if (!strcmp(*argv, "-num_general_map")) { 
        argc--; argv++; num_photons = atoi(*argv); 
      } else if (!strcmp(*argv, "-num_caustic_map")) { 
        argc--; argv++; num_caustics = atoi(*argv); 
      } else if (!strcmp(*argv, "-num_photon_estimate")) { 
        argc--; argv++; num_photon_estimate = atoi(*argv); 
      } else if (!strcmp(*argv, "-tone_map_const")) { 
        argc--; argv++; tone_map_const = atof(*argv); 
      } else { 
        fprintf(stderr, "Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      if (!input_scene_name) input_scene_name = *argv;
      else if (!output_image_name) output_image_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check scene filename
  if (!input_scene_name) {
    fprintf(stderr, "Usage: photonmap inputscenefile [outputimagefile] [-resolution <int> <int>] [-v]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// photons                                                            //
////////////////////////////////////////////////////////////////////////

// get's axis of rotation a theta to rotate vector a to b
void getConversionRotation(R3Vector a, R3Vector b, RNScalar *rotation_angle, R3Vector *rotation_axis)
{
    a.Normalize();
    RNScalar dot = a.Dot(b);
    R3Vector conv_rotation_axis = a;
    conv_rotation_axis.Cross(b);
    RNScalar conv_angle = acos(dot);
    if (RNIsNegative(conv_rotation_axis.Dot(b))) {
      conv_angle = RN_TWO_PI - conv_angle;
    }
    *rotation_angle = conv_angle;
    *rotation_axis = conv_rotation_axis;
}

void getR3CircleAxes(R3Vector normal, R3Vector *axis1, R3Vector *axis2)
{
  // Get circle axes
      RNDimension dim = normal.MinDimension();
      *axis1 = normal % R3xyz_triad[dim];
      axis1->Normalize();
      *axis2 = normal % *axis1;
      axis2->Normalize();
}
R3Point GetPhotonPosition(Photon *photon, void *dummy)
{
  return photon->position;
}

RNScalar getInteractionProbability(RNRgb power, RNRgb coeff)
{ 
  RNScalar numer = std::max({power[0] * coeff[0], power[1] * coeff[1], power[2] * coeff[2]});
  RNScalar denom = std::max({power[0], power[1], power[2]});
  return numer/denom;
}

void photonInteraction(const R3Brdf *brdf, RNScalar *prev_ior, const Photon* in, const R3Vector&  normal, RNRgb *out_power, R3Vector *out_direction, bool *absorbed, bool *is_diffuse, bool *is_specular_reflection)
{ 
  *absorbed = false;
  *is_diffuse = false;
  RNScalar diff = getInteractionProbability(in->power, brdf->Diffuse());
  RNScalar spec = getInteractionProbability(in->power, brdf->Specular());
  RNScalar trans = getInteractionProbability(in->power, brdf->Transmission());
  RNScalar absorb = 1 - diff - spec - trans;

  // scale probabilities if material does not objey conservation of 
  // energy
  if (RNIsGreater(diff + spec + trans, 1)) {
    RNScalar scaling =  RNScalar(1) / (diff + spec + trans);
    diff *= scaling;
    spec *= scaling;
    trans *= scaling;
    absorb = 1 - diff - spec - trans;
    assert(RNIsEqual(absorb, 0)); 
  } 
  assert(RNIsEqual(diff + spec + trans + absorb, 1));

  // Determine type of bounce and set new direction and power.
  RNScalar ksi = RNRandomScalar(); // random variable ksi
  if (ksi < diff) {
    // DIFFUSE case
    *out_power = (in->power * brdf->Diffuse()) / diff; 
    *is_diffuse = true;
    R3Vector rotation_axis;
    RNScalar rotation_angle;
    getConversionRotation(R3Vector(0,0,1), normal, &rotation_angle, &rotation_axis);

    RNScalar z = sqrt(RNRandomScalar());
    RNScalar  phi = RN_TWO_PI * RNRandomScalar();
    RNScalar theta = acos(z); 
    RNScalar sin_theta = sin(theta);
    RNScalar x = sin_theta * cos(phi);
    RNScalar y = sin_theta * sin(phi);
    R3Vector dir = R3Vector(x,y,z);
    dir.Rotate(rotation_axis, rotation_angle);
    dir.Normalize();
    *out_direction = dir;
  }
  else if (ksi <= diff + spec) {
    // SPECULAR case
    *out_power = (in->power * brdf->Specular()) / spec;
    *is_specular_reflection = true;
    R3Vector rotation_axis;
    RNScalar rotation_angle;
    R3Vector ideal_reflection = in->direction - (2 * in->direction.Dot(normal) * normal);
    ideal_reflection.Normalize();
    getConversionRotation(R3Vector(0,0,1), ideal_reflection, &rotation_angle, &rotation_axis);
    R3Vector dir;
    if (RNIsLess(ideal_reflection.Dot(normal),0)) {
      *absorbed = true;
      return;
    }
    do {
      RNScalar z = pow(RNRandomScalar(), RNScalar(1)/ brdf->Shininess());
      RNScalar alpha = acos(z); // angle between refleccted and ideal specular ray 
      RNScalar phi = RNRandomScalar() * RN_TWO_PI;
      RNScalar sin_alpha = sin(alpha);
      RNScalar x = sin_alpha * cos(phi);
      RNScalar y = sin_alpha * sin(phi);
      dir = R3Vector(x,y,z);
      dir.Rotate(rotation_axis, rotation_angle);
      dir.Normalize();
    } while (RNIsLess(dir.Dot(normal), 0));
    *out_direction = dir;
  } else if (ksi <= diff + spec + trans) {
  // TRANSMITTED case
    *out_power = (in->power * brdf->Transmission()) / trans;
    RNScalar ior_ratio = (brdf->IndexOfRefraction() == *prev_ior) ?  *prev_ior / camera_index_of_refraction : *prev_ior / brdf->IndexOfRefraction();
    RNScalar cos_theta = normal.Dot(-in->direction);
    if (cos_theta < 0) {
      cos_theta = -cos_theta;
      ior_ratio = RNScalar(1) / ior_ratio;
    }
    RNScalar sin_2_theta_t = pow(ior_ratio, 2) * 1 - pow(cos_theta, 2);

    if (RNIsGreaterOrEqual(sin_2_theta_t, RNScalar(1))) {
      // Total internal reflection
      *out_direction = in->direction - (2 * cos_theta * normal);
      return;
    }
    *prev_ior = brdf->IndexOfRefraction();
    *out_direction = (in->direction * ior_ratio) + ((ior_ratio * cos_theta - sqrt(1 - sin_2_theta_t)) * normal);
  }
   else {
    // ABSORBED case
    *absorbed = true;
  }
}

void tracePhoton(R3Scene *scene, RNScalar *prev_ior, Photon *in_photon,  RNArray<Photon *> &photon_list, bool is_caustic_map)
{
  // randomly terminate to prevent infinite photon tracing
  if (RNRandomScalar() < termination_rate) {
    return;
  }
  // Convenient variables
  R3SceneElement *element;
  R3Point point;
  R3Vector normal;

  in_photon->direction.Normalize();
  R3Ray ray = R3Ray(in_photon->source, in_photon->direction);

  bool is_bounce_allowed = in_photon->bounces < max_bounces;
  if (max_bounces == -1) {
    is_bounce_allowed = true;
  }

  if (!is_bounce_allowed || !(scene->Intersects(ray, NULL, &element, NULL, &point, &normal, NULL))) {
    return;
  }

  normal.Normalize();
  in_photon->normal = normal;
  in_photon->position = point;

  // Get intersection information
  const R3Material *material = (element) ? element->Material() : &R3default_material;
  const R3Brdf *brdf = (material) ? material->Brdf() : &R3default_brdf;

  RNRgb out_photon_power;
  R3Vector out_photon_direction;
  bool is_absorbed = false;
  bool is_diffuse = false;
  bool is_specular_reflection = false;
  photonInteraction(brdf, prev_ior, in_photon, normal, &out_photon_power, &out_photon_direction,  &is_absorbed, &is_diffuse, &is_specular_reflection);
  if (is_absorbed) {
    return;
  }
  Photon *out_photon = new Photon();

  if (is_diffuse && in_photon->bounces != 0) {
    photon_list.Insert(in_photon);
    if (is_caustic_map) {
      return;
    }
  }
  out_photon->direction = out_photon_direction;
  out_photon->direction.Normalize();
  // displace source slightly to avoid intersecting with same surface due to floating point error
  out_photon->source = point + RN_EPSILON * out_photon->direction;
  out_photon->position = point;
  out_photon->power = out_photon_power;
  out_photon->bounces = in_photon->bounces + 1;

  tracePhoton(scene, prev_ior, out_photon, photon_list, is_caustic_map);
}

// returns true if ray's first intersection is a specular reflection or transmission
bool IsRaySpecular(R3Scene *scene, R3Ray ray)
{
  R3Point point;
  R3Vector normal;
  R3SceneElement *element;
  if (!(scene->Intersects(ray, NULL, &element, NULL, &point, &normal, NULL))) {
    return false;
  }
  // Get intersection information
  const R3Material *material = (element) ? element->Material() : &R3default_material;
  const R3Brdf *brdf = (material) ? (material->Brdf()) : (&R3default_brdf);

  bool is_absorbed = false;
  bool is_diffuse = false;
  bool is_specular_reflection = false;

  Photon in_photon;
  in_photon.direction = ray.Vector();
  in_photon.direction.Normalize();
  in_photon.power = RNRgb(1,1,1);
  RNRgb out_photon_power;
  R3Vector out_photon_direction;
  RNScalar prev_ior = camera_index_of_refraction;
  photonInteraction(brdf, &prev_ior, &in_photon, normal, &out_photon_power, &out_photon_direction,  &is_absorbed, &is_diffuse, &is_specular_reflection);
  if (is_diffuse || is_absorbed) {
    return false;
  }
  return true;
}

// reccursively traces ray until a diffuse interaction with a surface
bool traceRayDiffuse(R3Scene *scene, RNScalar *prev_ior, R3Ray ray, R3Point *point, R3SceneElement **element , R3Vector *normal, RNScalar termination_rate_ray_trace, RNRgb *power_multiplier)
{ 
  // randomly terminate to prevent infinite photon tracing
  if (RNRandomScalar() < termination_rate_ray_trace) {
    return false;
  }

  if (!(scene->Intersects(ray, NULL, element, NULL, point, normal, NULL))) {
    return false;
  }
  normal->Normalize();
  // Get intersection information
  const R3Material *material = (*element) ? (*element)->Material() : &R3default_material;
  const R3Brdf *brdf = (material) ? (material->Brdf()) : (&R3default_brdf);

  bool is_absorbed = false;
  bool is_diffuse = false;
  bool is_specular_reflection = false;

  Photon in_photon;
  in_photon.direction = ray.Vector();
  in_photon.direction.Normalize();
  in_photon.power = RNRgb(1,1,1);
  RNRgb out_photon_power;
  R3Vector out_photon_direction;
  photonInteraction(brdf, prev_ior, &in_photon, *normal, &out_photon_power, &out_photon_direction,  &is_absorbed, &is_diffuse, &is_specular_reflection);
  out_photon_direction.Normalize();
  if (is_absorbed) {
    return false;
  } else if (is_diffuse) {
    *power_multiplier = *power_multiplier * out_photon_power;
    return true;
  } else if (is_specular_reflection) {
    // adjust for importance sampling
    *power_multiplier = *power_multiplier * out_photon_power * (brdf->Shininess() +2) / (brdf->Shininess() + 1);
  }
  ray = R3Ray(*point + RN_EPSILON * out_photon_direction, out_photon_direction, false);
  return traceRayDiffuse(scene, prev_ior, ray, point, element, normal, termination_rate_ray_trace, power_multiplier);
}

RNRgb EstimateFlux(R3Kdtree<Photon *> *photon_map,  R3Point point, int num_photons, RNScalar max_distance, RNRgb diffuseBrdf) {
  RNArray<Photon *> nearby;
  RNLength *distances = new RNLength[num_photons];
  photon_map->FindClosest(point, RNScalar(0), max_distance, num_photons, nearby, distances);
  if (nearby.NEntries() == 0) {
    return RNRgb(0,0,0);
  }
  RNRgb color = RNRgb(0,0,0);
  int num_nearby = nearby.NEntries();
  for (int i = 0; i < num_nearby; ++i)
  {
    color += nearby[i]->power;
  }
  color *= diffuseBrdf;
  color /= RN_PI * distances[num_nearby - 1] * distances[num_nearby - 1];
  delete [] distances;
  return color;
}

// initilize photons on light source
static RNArray<Photon *>  
GetPhotonsFromLights(R3Scene *scene, long num_photons, bool get_only_caustics)
{
  RNArray<Photon *>  photons_from_lights;
  // num_photons is total number of photons emitted from the lights that 
  // intersect with the secne.


  int total_intensity = 0;
  for (int k = 0; k < scene->NLights(); k++) {
    R3Light *light = scene->Light(k);
    total_intensity += light->Intensity();
  }

  long photons_per_intesity = (num_photons + num_caustics) / total_intensity;
  RNScalar photon_power = RNScalar(1)/photons_per_intesity;

  // Draw intersection point and normal for some rays
  for (int k = 0; k < scene->NLights(); k++) {
    R3Light *light = scene->Light(k);

     if (light->ClassID() == R3PointLight::CLASS_ID()) {
      // Point light case
      R3PointLight *point_light = (R3PointLight *) light;
      for (long i = 0; i<photons_per_intesity * point_light->Intensity(); i++) {
        Photon *curr_photon = new Photon();
        R3Ray ray;
        do {
          double x = 2.0 * RNRandomScalar() - 1.0;
          double y = 2.0 * RNRandomScalar() - 1.0;
          double z = 2.0 * RNRandomScalar() - 1.0;
          curr_photon->direction = R3Vector(x, y, z);
          ray = R3Ray(point_light->Position(), curr_photon->direction);
          curr_photon->normal = R3Vector(0,0,0);
          curr_photon->source = point_light->Position();
          curr_photon->position = R3Point(0,0,0); // to initilize position
          curr_photon->bounces = 0;
          curr_photon->power  = point_light->Color() * photon_power;
        } while (curr_photon->direction.Length() > 1);
        if (!get_only_caustics || (get_only_caustics && IsRaySpecular(scene, R3Ray(curr_photon->source, curr_photon->direction))))
        photons_from_lights.Insert(curr_photon);
      }
    } else if (light->ClassID() == R3SpotLight::CLASS_ID()) {
      // Point light case
      R3SpotLight *spot_light = (R3SpotLight *) light;
      for (long i = 0; i<photons_per_intesity * spot_light->Intensity(); i++) {
        Photon *curr_photon = new Photon();
        R3Ray ray;

        // sample points in the cone defined by cone
        R3Point position = spot_light->Position();
        R3Vector central_direction = spot_light->Direction();
        R3Vector rotation_axis;
        RNScalar rotation_angle;
        central_direction.Normalize();
        getConversionRotation(R3Vector(0,0,1), central_direction, &rotation_angle, &rotation_axis);
        RNScalar bias_towards_center = 3;
        do {
          RNScalar z = pow(RNRandomScalar(), RNScalar(1)/ bias_towards_center);
          RNScalar alpha = acos(z); // angle between refleccted and ideal specular ray 
          RNScalar phi = RNRandomScalar() * RN_TWO_PI;
          RNScalar sin_alpha = sin(alpha);
          RNScalar x = sin_alpha * cos(phi);
          RNScalar y = sin_alpha * sin(phi);
          R3Vector dir = R3Vector(x,y,z);
          dir.Rotate(rotation_axis, rotation_angle);
          dir.Normalize();

          curr_photon->direction = dir;
          ray = R3Ray(spot_light->Position(), curr_photon->direction);
          curr_photon->normal = R3Vector(0,0,0);
          curr_photon->source = spot_light->Position();
          curr_photon->position = R3Point(0,0,0); // to initilize position
          curr_photon->bounces = 0;
          curr_photon->power = spot_light->Color() * photon_power;
        } while (curr_photon->direction.Dot(central_direction) < cos(spot_light->CutOffAngle()));
        if (!get_only_caustics || (get_only_caustics && IsRaySpecular(scene, R3Ray(curr_photon->source, curr_photon->direction))))
        photons_from_lights.Insert(curr_photon);
      }
    } else if (light->ClassID() == R3DirectionalLight::CLASS_ID()) {
      // directional light case
      double radius = scene->BBox().DiagonalRadius();
      R3Point scene_center = scene->BBox().Centroid();
      R3DirectionalLight *dir_light = (R3DirectionalLight *) light;
      R3Vector dir_light_dir = dir_light->Direction();
      dir_light_dir.Normalize();
      R3Point dir_light_pos = scene_center - (dir_light_dir * radius);
      R3Vector axis1;
      R3Vector axis2;
      getR3CircleAxes(dir_light_dir, &axis1, &axis2);
        
      for (long i = 0; i<photons_per_intesity * dir_light->Intensity(); i++) {
        Photon *curr_photon = new Photon();
        R3Ray ray;
        // sample points  on circle to shoot photons from
        R3Point source_pos;
        RNScalar r1;
        RNScalar r2;
        do {
          r1 = (RNRandomScalar() * 2) - 1;
          r2 = (RNRandomScalar() * 2) - 1;
        } while(r1 * r1 + r2 * r2 > 1);
        source_pos = dir_light_pos;
        source_pos += (r1 * axis1 * radius) + (r2 * axis2 * radius);
        curr_photon->direction = dir_light_dir;
        ray = R3Ray(source_pos, curr_photon->direction);
        curr_photon->normal = R3Vector(0,0,0);
        curr_photon->source = source_pos;
        curr_photon->position = R3Point(0,0,0);
        curr_photon->bounces = 0;
        curr_photon->power = dir_light->Color() * photon_power;
        if (!get_only_caustics || (get_only_caustics && IsRaySpecular(scene, R3Ray(curr_photon->source, curr_photon->direction))))
        photons_from_lights.Insert(curr_photon);
      }
    } else if (light->ClassID() == R3AreaLight::CLASS_ID()) {
      // Point light case
      R3AreaLight *area_light = (R3AreaLight *) light;
      R3Vector axis1;
      R3Vector axis2;
      getR3CircleAxes(area_light->Direction(), &axis1, &axis2);
      for (long i = 0; i<photons_per_intesity * area_light->Intensity(); i++) {
        Photon *curr_photon = new Photon();
        R3Ray ray;
        // sample points  on circle to shoot photons from
        R3Point source_pos;
        RNScalar r1;
        RNScalar r2;
        do {
          r1 = (RNRandomScalar() * 2) - 1;
          r2 = (RNRandomScalar() * 2) - 1;
        } while(r1 * r1 + r2 * r2 > 1);
        source_pos = area_light->Position();
        source_pos += (r1 * axis1 * area_light->Radius()) + (r2 * axis2 * area_light->Radius());

        R3Vector rotation_axis;
        RNScalar rotation_angle;
        getConversionRotation(R3Vector(0,0,1), area_light->Direction(), &rotation_angle, &rotation_axis);
        RNScalar z = sqrt(RNRandomScalar());
        RNScalar phi = RN_TWO_PI * RNRandomScalar();
        RNScalar theta = acos(z); 
        RNScalar sin_theta = sin(theta);
        RNScalar x = sin_theta * cos(phi);
        RNScalar y = sin_theta * sin(phi);
        R3Vector dir = R3Vector(x,y,z);
        dir.Rotate(rotation_axis, rotation_angle);
        dir.Normalize();

        curr_photon->direction = dir;
        ray = R3Ray(source_pos, curr_photon->direction);
        curr_photon->normal = R3Vector(0,0,0);
        curr_photon->source = source_pos;
        curr_photon->position = R3Point(0,0,0); 
        curr_photon->bounces = 0;
        curr_photon->power = area_light->Color() * photon_power;
        if (!get_only_caustics || (get_only_caustics && IsRaySpecular(scene, R3Ray(curr_photon->source, curr_photon->direction))))
        photons_from_lights.Insert(curr_photon);
      }
    } else {
      std::cout << "unrecognized light" << std::endl;
      assert(false);
      return photons_from_lights;
    } 
  }
  return photons_from_lights;
}
////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);



  // Read scene
  std::cout << input_scene_name << std::endl;
  scene = ReadScene(input_scene_name);
  if (!scene) exit(-1);

  // Initialize photons
  RNArray<Photon *> photons_from_lights = GetPhotonsFromLights(scene, num_photons, false); //  get_only_caustics  = false
  RNArray<Photon *> caustics_from_lights = GetPhotonsFromLights(scene, num_caustics, true); //  get_only_caustics = true

  std::cout<<"shooting photons..."<< std::endl;
  RNScalar russian_roulette_multiplier = RNScalar(1) / 1 - termination_rate;
  for (int i = 0; i < photons_from_lights.NEntries(); i++) {
    // N.B we assume that camera is in vaccum
    Photon *photon_from_light = photons_from_lights[i];
    photon_from_light->power *= russian_roulette_multiplier;
    RNScalar ior = camera_index_of_refraction;
    tracePhoton(scene, &ior, photon_from_light, photon_list, false); // is_caustic_map = false
  } 

  for (int i = 0; i < caustics_from_lights.NEntries(); i++) {
    // N.B we assume that camera is in vaccum
    Photon *caustic_from_light = caustics_from_lights[i];
    caustic_from_light->power *= russian_roulette_multiplier;
    RNScalar ior = camera_index_of_refraction;
    tracePhoton(scene, &ior, caustic_from_light, caustic_list, true); // is_caustic_map = true
  } 

  photon_map = new R3Kdtree<Photon *>(photon_list, GetPhotonPosition);
  caustic_map = new R3Kdtree<Photon *>(caustic_list, GetPhotonPosition);
  std::cout<<"photon mapping done, now rendering.."<< std::endl;
  // Check output image file
  if (output_image_name) {
    // Set scene viewport
    scene->SetViewport(R2Viewport(0, 0, render_image_width, render_image_height));
    // Render image
    R2Image *image = RenderImage(scene, photon_map, caustic_map, render_image_width, render_image_height, print_verbose, num_samples, general_search_range,caustic_search_range, tone_map_const, num_photon_estimate);
    if (!image) exit(-1);

    // Write image
    if (!WriteImage(image, output_image_name)) exit(-1);

    // Delete image
    delete image;
  }
  else {
    // Initialize GLUT
    GLUTInit(&argc, argv);

    // Create viewer
    viewer = new R3Viewer(scene->Viewer());
    if (!viewer) exit(-1);
    
    // Run GLUT interface
    GLUTMainLoop();
  }
  // Return success 
  return 0;
  // Delete viewer (doesn't ever get here)
  delete viewer;
  delete photon_map;
  for (int i = 0; i < photon_list.NEntries(); i++) {
    delete photon_list[i];
  }
  delete caustic_map;
  for (int i = 0; i < caustic_list.NEntries(); i++) {
    delete caustic_list[i];
  }
}

















