/*
  System includes
*/

#include <stdio.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>

/*
  Debug level:
  Comment out or set to 0 for no output
  Set to 1 for everything but suspendcheck
  Set to 2 for everything
*/

#define DEBUG 1

/*
  Function prototypes
*/

typedef void (*sighandler_t)(int);
sighandler_t signal(int signum, sighandler_t handler);
static void suspendlog(char *msg);
static sighandler_t suspendhandler();

/*
  Local variables
*/

static int needSuspend=0;

/*
  Logs messages to stdout if DEBUG is set
*/

static void inline suspendlog(char *msg)
{
#if DEBUG
  printf("%s ; pid=%d;needSuspend=%d\n",msg,getpid(),needSuspend); 
  fflush(stdout);
#endif
}

/*
  Signal handler.
  Called if we get a TSTP signal.
  Just sets needSuspend.
  suspendcheck will suspend our process next time called after needSuspend is set 
*/
static sighandler_t suspendhandler(int ignore)
{
  needSuspend=1;
  suspendlog("suspendhandler called");
}
/*
  Initialization routine.
  Should be called once by user, after MPIinitialize
  Activates the TSTP signal handler
*/
void suspendinit_() {
  suspendlog("suspendinit called");
  signal(SIGTSTP,(sighandler_t) suspendhandler);
}
/*
  Suspend check
  Should be called regularly by user, in 'safe', compute-only code sections.
*/
void suspendcheck_() {
#if DEBUG > 1
  suspendlog("suspendcheck called");
#endif
  if(needSuspend) {
    kill(getpid(),SIGSTOP);
    needSuspend=0;
  }
}
