#include <signal.h>
#include "fem.h"

CWManager *CWManager::instance = 0;


WINDOW *create_newwin(int height, int width, int starty, int startx)
{
  WINDOW *local_win;

  local_win = newwin(height, width, starty, startx);
  box(local_win, 0 , 0);  /* 0, 0 gives default characters 
			   * for the vertical and horizontal
			   * lines			    */
  wrefresh(local_win);    /* Show that box 	            */
  
  delwin(local_win);
  local_win = newwin(height-2, width-2, starty+1, startx+1);
						
  return local_win;
}
						    
void destroy_win(WINDOW *local_win)
{
  /* box(local_win, ' ', ' '); : This won't produce the desired
   * result of erasing the window. It will leave it's four corners 
   * and so an ugly remnant of window. 
   */
  wborder(local_win, ' ', ' ', ' ',' ',' ',' ',' ',' ');
  /* The parameters taken are 
   * 1. win: the window on which to operate
   * 2. ls: character to be used for the left side of the window 
   * 3. rs: character to be used for the right side of the window 
   * 4. ts: character to be used for the top side of the window 
   * 5. bs: character to be used for the bottom side of the window 
   * 6. tl: character to be used for the top left corner of the window 
   * 7. tr: character to be used for the top right corner of the window 
   * 8. bl: character to be used for the bottom left corner of the window 
   * 9. br: character to be used for the bottom right corner of the window
   */
  wrefresh(local_win);
  delwin(local_win);
}

void terminate(int s)
{
  if (CWManager::getInstanceNotCreate() != 0)
  {
    wprintw(CWManager::getInstance()->getProgW(), "\nProgram terminated.\n");
    delete CWManager::getInstance();
  }
  exit(s);
}


CWManager *CWManager::getInstance()
{
  if (instance == 0)
  {
    instance = new CWManager;
  }
  
  return instance;
}

CWManager *CWManager::getInstanceNotCreate()
{
  return instance;
}

CWManager::CWManager()
{
  // initialize ncurses
  initscr();
  
  // get immediately pressed keys, signals (Ctrl+C, Ctrl+Z) are working.
  cbreak();    
  
  // handle signals
  (void)signal(SIGINT, terminate);
  
  // suppress printing of typed text
  noecho();

  // print program title
  attron(A_BOLD);
  mvprintw(0, 0, "%s v. %s", PROGRAM_NAME, PROGRAM_VERSION);
  attroff(A_BOLD);
  refresh();

  // create main window
  progWx = 0;
  progWy = 1;
  progWw = COLS/2;
  progWh = LINES - progWy;
  progW = create_newwin(progWh, progWw, progWy, progWx);
  
  // create solver window
  solveWx = progWw;
  solveWy = progWy;
  solveWw = COLS - progWw;
  solveWh = LINES - solveWy;
  solveW = create_newwin(solveWh, solveWw, solveWy, solveWx);
  
  // enable window scrolling
  scrollok(progW, true);
  scrollok(solveW, true);
}



CWManager::~CWManager()
{
  wattron(progW, A_BOLD);
  wprintw(progW, "\nPress any key to exit.");
  wattroff(progW, A_BOLD);
  wscanw(progW, "");

  destroy_win(progW);
  destroy_win(solveW);
  
  endwin();
  
  instance = 0;
}

