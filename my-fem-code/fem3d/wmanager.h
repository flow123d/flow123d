#ifndef _WMANAGER_H__
#define _WMANAGER_H__

#include <ncurses.h>

#define W_PROG  1
#define W_SOLVE 2

WINDOW *create_newwin(int height, int width, int starty, int startx);
void    destroy_win(WINDOW *local_win);

/** Singleton class for handling windows by ncurses. */
class CWManager
{
 private:
  static CWManager *instance; /** Instance of window manager. */
 
  WINDOW *progW,  /** Window with main program information. */
         *solveW; /** Window with solution information. */
  int progWx, progWy, progWw, progWh,
      solveWx, solveWy, solveWw, solveWh; /**< Window dimensions. */
  
 public:
  CWManager();
  ~CWManager();
  
  /** Returns pointer to WM (and creates it if it does not exist). */
  static CWManager *getInstance();
  
  /** Returns pointer to WM whether it exists or not. */
  static CWManager *getInstanceNotCreate();
  
  WINDOW *getProgW() { return progW; };
  WINDOW *getSolveW() { return solveW; };
};


#endif
