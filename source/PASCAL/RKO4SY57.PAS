program ALG057;
{     RUNGE-KUTTA METHOD FOR SYSTEMS OF DIFFERENTIAL EQUATIONS

      TO APPROXIMATE THE SOLUTION OF THE MTH-ORDER SYSTEM OF FIRST-
      ORDER INITIAL-VALUE PROBLEMS
                 UJ' = FJ( T, U1, U2, ..., UM ), J = 1, 2, ..., M
                  A <= T <= B, UJ(A) = ALPHAJ, J = 1, 2, ..., M
      AT (N+1) EQUALLY SPACED NUMBERS IN THE INTERVAL A <= T <= B.

      INPUT:   ENDPOINTS A,B; NUMBER OF EQUATIONS M; INTEGER N; INITIAL
               CONDITIONS ALPHA1, ..., ALPHAM.

      OUTPUT:  APPROXIMATION WJ TO UJ(T) AT THE (N+1) VALUES OF T.
}
var
   A,B,ALPHA1,ALPHA2,H,T,W1,W2,X11,X12,X21,X22,X31,X32,X41,X42 : real;
   FLAG,N,I : integer;
   AA : char;
   OK : boolean;
   NAME : string [ 14 ];
   OUP : text;
function F1 ( T, X1, X2 : real ) : real;
   begin
      F1 := -4*X1+3*X2+6
   end;
function F2 ( T, X1, X2 : real ) : real;
   begin
      F2 := -2.4*X1+1.6*X2+3.6
   end;
procedure INPUT;
begin
   writeln('This is the Runge-Kutta Method for Systems with m = 2.');
   OK := false;
   write ('Have the functions F1 and F2 been defined? ');
   writeln ('Answer Y or N. ');
   readln ( AA );
   if ( AA = 'Y' ) or ( AA = 'y' ) then
      begin
         OK := false;
         while ( not OK ) do
            begin
               write ('Input left and right endpoints separated by blank ');
               writeln ('- include decimal point ');
               readln ( A, B );
               if ( A >= B ) then
                  writeln ('Left endpoint must be less than right endpoint ')
               else OK := true
            end;
         OK := false;
         while ( not OK ) do
            begin
               write ('Input a positive integer for the number of ');
               writeln ('subintervals ');
               readln ( N );
               if ( N <= 0 ) then
                  writeln ('Number must be a positive integer ')
               else OK := true
            end;
         write ('Input the two initial conditions, seperated by blank ');
         writeln ('- include decimal point ');
         readln ( ALPHA1, ALPHA2 )
      end
   else
      writeln ('The program will end so that the functions can be created.')
end;
procedure OUTPUT;
   begin
      writeln ('Choice of output method: ');
      writeln ('1. Output to screen ');
      writeln ('2. Output to text file ');
      writeln ('Please enter 1 or 2 ');
      readln ( FLAG );
      if ( FLAG = 2 ) then
         begin
            writeln ('Input the file name in the form - drive:name.ext ');
            readln ( NAME );
            assign ( OUP, NAME )
         end
      else assign ( OUP, 'CON' );
      rewrite ( OUP );
      writeln(OUP,'RUNGE-KUTTA METHOD FOR SYSTEMS WITH m = 2.');
      writeln ( OUP, 'T':5,'W1':12,'W2':12);
      writeln ( OUP );
   end;
begin
   INPUT;
   if OK then
      begin
         OUTPUT;
{     STEP 1                                                                   }
         H := ( B - A ) / N;
         T := A;
{     STEP 2                                                                   }
         W1 := ALPHA1;
         W2 := ALPHA2;
{     STEP 3                                                                   }
         writeln ( OUP,T:5:3,W1:12:8,W2:12:8);
{     STEP 4                                                                   }
         for I := 1 to N do
            begin
{           STEP 5                                                             }
               X11 := H * F1( T, W1, W2 );
               X12 := H * F2( T, W1, W2 );
{           STEP 6                                                             }
               X21 := H * F1( T + H / 2.0, W1 + X11 / 2.0, W2 + X12 / 2.0 );
               X22 := H * F2( T + H / 2.0, W1 + X11 / 2.0, W2 + X12 / 2.0 );
{           STEP 7                                                             }
               X31 := H * F1( T + H / 2.0, W1 + X21 / 2.0, W2 + X22 / 2.0 );
               X32 := H * F2( T + H / 2.0, W1 + X21 / 2.0, W2 + X22 / 2.0 );
{           SYEP 8                                                             }
               X41 := H * F1( T + H, W1 + X31, W2 + X32 );
               X42 := H * F2( T + H, W1 + X31, W2 + X32 );
{           STEP 9                                                             }
               W1 := W1 + ( X11 + 2.0 * X21 + 2.0 * X31 + X41 ) / 6.0;
               W2 := W2 + ( X12 + 2.0 * X22 + 2.0 * X32 + X42 ) / 6.0;
{           STEP 10                                                            }
               T := A + I * H;
{           STEP 11                                                            }
               writeln ( OUP,T:5:3,W1:12:8,W2:12:8);
            end;
{     STEP 12                                                                  }
         close ( OUP )
      end
end.
