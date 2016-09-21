unit LUX.DiscreteTrans.D1;

interface //#################################################################### Å°

uses LUX;

type //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Åyå^Åz

     TDiscreteTrans = class;

     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyÉåÉRÅ[ÉhÅz

     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyÉNÉâÉXÅz

     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TDiscreteTrans

     TDiscreteTrans = class
     private
       _NC :Integer;
       ///// initializing routines
       procedure makeipt;
       procedure makewt;
       procedure makect( var c:array of Double );
       /////
       procedure bitrv2;
       procedure bitrv2conj;

       procedure bitrv216;
       procedure bitrv216neg;

       procedure bitrv208;
       procedure bitrv208neg;

       procedure cftf1st( var w:array of Double );
       procedure cftb1st( var w:array of Double );

       procedure cftmdl1( const n:Integer; var a:array of Double; var w:array of Double );
       procedure cftmdl2( const n:Integer; var a:array of Double; var w:array of Double );

       function cfttree( const n,j,k:Integer ) :Integer;

       procedure cftf161( var a:array of Double; var w:array of Double );
       procedure cftf162( var a:array of Double; var w:array of Double );

       procedure cftf081( var a:array of Double; var w:array of Double );
       procedure cftf082( var a:array of Double; var w:array of Double );

       procedure cftleaf( const n,isplt:Integer; var a:array of Double );

       procedure cftrec4;

       procedure cftfx41;

       procedure cftf040;
       procedure cftb040;

       procedure cftx020;
     protected
       _TEMP  :array of Double;
       _TempN :Integer;
       _IP    :array of Integer;
       _W     :array of Double;
       _NW    :Integer;
       _NormW :Double;
       ///// ÉtÉBÅ[ÉãÉh
       _Count :Integer;  upCount :Boolean;
       ///// ÉAÉNÉZÉX
       procedure SetCount( const Count_:Integer ); virtual;
       ///// ÉÅÉ\ÉbÉh
       procedure cftfsub;
       procedure cftbsub;
       procedure rftfsub( var c:array of Double );
       procedure rftbsub( var c:array of Double );

       procedure Normalize;

       procedure MakeTableW;
       procedure MakeTableC;
     public
       constructor Create;
       destructor Destroy; override;
       ///// ÉvÉçÉpÉeÉB
       property Count :Integer read _Count write SetCount;
       ///// ÉÅÉ\ÉbÉh
       procedure TransWF; virtual; abstract;
       procedure TransFW; virtual; abstract;
     end;

//const //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyíËêîÅz

//var //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyïœêîÅz

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyÉãÅ[É`ÉìÅz

implementation //################################################################################### Å°

uses Math;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyÉåÉRÅ[ÉhÅz

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyÉNÉâÉXÅz

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TDiscreteTrans

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& private

////////////////////////////////////////////////////////////////////////////// initializing routines

procedure TDiscreteTrans.makeipt;
var
   j, l, m, m2, p, q :Integer;
begin
     _IP[ 2 ] :=  0;
     _IP[ 3 ] := 16;

     m := 2;
     l := _NW;
     while l > 32 do
     begin
          m2 := m  shl 1;
          q  := m2 shl 3;

          for j := m to m2 - 1 do
          begin
               p := _IP[ j ] shl 2;

               _IP[ m  + j ] := p;
               _IP[ m2 + j ] := p + q;
          end;

          m := m2;
          l := l shr 2
     end
end;

procedure TDiscreteTrans.makewt;
var
   j, nwh, nw0, nw1 :Integer;
   delta, wn4r, wk1r, wk1i, wk3r, wk3i :Double;
begin
     _IP[ 0 ] := _NW;
     _IP[ 1 ] := 1;

     if _NW > 2 then
     begin
          nwh   := _NW shr 1;
          delta := ArcTan( 1 ) / nwh;
          wn4r  := Cos( delta * nwh );

          _W[ 0 ] := 1;
          _W[ 1 ] := wn4r;

          if nwh = 4 then
          begin
               _W[ 2 ] := Cos( delta * 2 );
               _W[ 3 ] := Sin( delta * 2 );
          end
          else
          if nwh > 4 then
          begin
               makeipt;

               _W[ 2 ] := 0.5 / Cos( delta * 2 );
               _W[ 3 ] := 0.5 / Cos( delta * 6 );

               j := 4;
               while j < nwh do
               begin
                    _W[ j     ] := +Cos(     delta * j );
                    _W[ j + 1 ] := +Sin(     delta * j );
                    _W[ j + 2 ] := +Cos( 3 * delta * j );
                    _W[ j + 3 ] := -Sin( 3 * delta * j );

                    Inc( j, 4 )
               end
          end;

          nw0 := 0;
          while nwh > 2 do
          begin
               nw1 := nw0 + nwh;
               nwh := nwh shr 1;

               _W[ nw1     ] := 1;
               _W[ nw1 + 1 ] := wn4r;

               if nwh = 4 then
               begin
                    wk1r := _W[ nw0 + 4 ];
                    wk1i := _W[ nw0 + 5 ];

                    _W[ nw1 + 2 ] := wk1r;
                    _W[ nw1 + 3 ] := wk1i;
               end
               else
               if nwh > 4 then
               begin
                    wk1r := _W[ nw0 + 4 ];
                    wk3r := _W[ nw0 + 6 ];

                    _W[ nw1 + 2 ] := 0.5 / wk1r;
                    _W[ nw1 + 3 ] := 0.5 / wk3r;

                    j := 4;
                    while j < nwh do
                    begin
                         wk1r := _W[ nw0 + 2 * j     ];
                         wk1i := _W[ nw0 + 2 * j + 1 ];
                         wk3r := _W[ nw0 + 2 * j + 2 ];
                         wk3i := _W[ nw0 + 2 * j + 3 ];

                         _W[ nw1 + j     ] := wk1r;
                         _W[ nw1 + j + 1 ] := wk1i;
                         _W[ nw1 + j + 2 ] := wk3r;
                         _W[ nw1 + j + 3 ] := wk3i;

                         Inc( j, 4 )
                    end
               end;

               nw0 := nw1
          end
     end
end;

procedure TDiscreteTrans.makect( var c:array of Double );
var
   j, nch :Integer;
   delta :Double;
begin
     _IP[ 1 ] := _NC;

     if _NC > 1 then
     begin
          nch   := _NC shr 1;
          delta := ArcTan( 1 ) / nch;

          c[  0  ] := Cos( delta * nch );
          c[ nch ] := 0.5 * c[ 0 ];

          for j := 1 to nch - 1do
          begin
               c[       j ] := 0.5 * Cos( delta * j );
               c[ _NC - j ] := 0.5 * Sin( delta * j );
          end
     end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.bitrv2;
var
   j, j1, k, k1, l, m, nh, nm :Integer;
   xr, xi, yr, yi :Double;
begin
     m := 1;
     l := _TempN shr 2;
     while l > 8 do
     begin
          m := m shl 1;
          l := l shr 2
     end;

     nh := _TempN shr 1;
     nm := 4 * m;
     if l = 8 then
     begin
          for k := 0 to m - 1 do
          begin
               for j := 0 to k - 1 do
               begin
                    j1 := 4 * j + 2 * _IP[ m + k ];
                    k1 := 4 * k + 2 * _IP[ m + j ];

                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Dec( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nh );
                    Inc( k1,  2 );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1,     nm );
                    Dec( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Inc( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1,     nm );
                    Dec( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,  2 );
                    Inc( k1, nh );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Dec( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nh );
                    Dec( k1,  2 );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Dec( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Inc( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Dec( k1, 2 * nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;
               end;

               k1 := 4 * k + 2 * _IP[ m + k ];
               j1 := k1 + 2;
               Inc( k1, nh );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

               Inc( j1, nm );
               Inc( k1, 2 * nm );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

               Inc( j1, nm );
               Dec( k1, nm );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

               Dec( j1, 2 );
               Dec( k1, nh );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

               Inc( j1, nh + 2 );
               Inc( k1, nh + 2 );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

               Dec( j1, nh - nm );
               Inc( k1, 2 * nm - 2 );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;
          end
     end
     else
     begin
          for k := 0 to m - 1 do
          begin
               for j := 0 to k - 1 do
               begin
                    j1 := 4 * j + _IP[ m + k ];
                    k1 := 4 * k + _IP[ m + j ];
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Inc( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nh );
                    Inc( k1, 2 );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Dec( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, 2 );
                    Inc( k1, nh );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Inc( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nh );
                    Dec( k1, 2 );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Dec( k1, nm );
                    xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
                    yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;
               end;

               k1 := 4 * k + _IP[ m + k ];
               j1 := k1 + 2;
               Inc( k1, nh );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;

               Inc( j1, nm );
               Inc( k1, nm );
               xr := _TEMP[ j1 ]; xi := _TEMP[ j1 + 1 ];
               yr := _TEMP[ k1 ]; yi := _TEMP[ k1 + 1 ];
               _TEMP[ j1 ] := yr; _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1 ] := xr; _TEMP[ k1 + 1 ] := xi;
          end
     end
end;

procedure TDiscreteTrans.bitrv2conj;
var
   j, j1, k, k1, l, m, nh, nm :Integer;
   xr, xi, yr, yi :Double;
begin
     m := 1;
     l := _TempN shr 2;
     while l > 8 do
     begin
          m := m shl 1;
          l := l shr 2;
     end;

     nh := _TempN shr 1;
     nm := 4 * m;
     if l = 8 then
     begin
          for k := 0 to m - 1 do
          begin
               for j := 0 to k - 1 do
               begin
                    j1 := 4 * j + 2 * _IP[ m + k ];
                    k1 := 4 * k + 2 * _IP[ m + j ];

                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1 ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1 ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Dec( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nh );
                    Inc( k1,  2 );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1,     nm );
                    Dec( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Inc( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1,     nm );
                    Dec( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,  2 );
                    Inc( k1, nh );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Dec( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,     nm );
                    Inc( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nh );
                    Dec( k1,  2 );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1,     nm );
                    Dec( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Inc( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1,     nm );
                    Dec( k1, 2 * nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;
               end;

               k1 := 4 * k + 2 * _IP[ m + k ];

               j1 := k1 + 2;
               Inc( k1, nh );
               _TEMP[ j1 - 1 ] := -_TEMP[ j1 - 1 ];
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;
               _TEMP[ k1 + 3 ] := -_TEMP[ k1 + 3 ];

               Inc( j1, nm );
               Inc( k1, 2 * nm );
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;

               Inc( j1, nm );
               Dec( k1, nm );
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;

               Dec( j1, 2 );
               Dec( k1, nh );
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;

               Inc( j1, nh + 2 );
               Inc( k1, nh + 2 );
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1    ] := xr;
               _TEMP[ k1 + 1 ] := xi;

               Dec( j1, nh - nm );
               Inc( k1, 2 * nm - 2 );
               _TEMP[ j1 - 1 ] := -_TEMP[ j1 - 1 ];
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;
               _TEMP[ k1 + 3 ] := -_TEMP[ k1 + 3 ];
          end
     end
     else
     begin
          for k := 0 to m - 1 do
          begin
               for j := 0 to k - 1 do
               begin
                    j1 := 4 * j + _IP[ m + k ];
                    k1 := 4 * k + _IP[ m + j ];

                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Inc( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nh );
                    Inc( k1,  2 );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Dec( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1,  2 );
                    Inc( k1, nh );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Inc( j1, nm );
                    Inc( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nh );
                    Dec( k1,  2 );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;

                    Dec( j1, nm );
                    Dec( k1, nm );
                    xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
                    yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
                    _TEMP[ j1     ] := yr;
                    _TEMP[ j1 + 1 ] := yi;
                    _TEMP[ k1     ] := xr;
                    _TEMP[ k1 + 1 ] := xi;
               end;

               k1 := 4 * k + _IP[ m + k ];
               j1 := k1 + 2;
               Inc( k1, nh );
               _TEMP[ j1 - 1 ] := -_TEMP[ j1 - 1 ];
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;
               _TEMP[ k1 + 3 ] := -_TEMP[ k1 + 3 ];

               Inc( j1, nm );
               Inc( k1, nm );
               _TEMP[ j1 - 1 ] := -_TEMP[ j1 - 1 ];
               xr := +_TEMP[ j1     ]; xi := -_TEMP[ j1 + 1 ];
               yr := +_TEMP[ k1     ]; yi := -_TEMP[ k1 + 1 ];
               _TEMP[ j1     ] := yr;
               _TEMP[ j1 + 1 ] := yi;
               _TEMP[ k1     ] := xr;
               _TEMP[ k1 + 1 ] := xi;
               _TEMP[ k1 + 3 ] := -_TEMP[ k1 + 3 ];
          end
     end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.bitrv216;
var
   x01r, x01i,
   x02r, x02i,
   x03r, x03i,
   x04r, x04i,
   x05r, x05i,
   {--}  {--}
   x07r, x07i,
   x08r, x08i,
   {--}  {--}
   x10r, x10i,
   x11r, x11i,
   x12r, x12i,
   x13r, x13i,
   x14r, x14i
   {--}  {--} :Double;
begin
    x01r := _TEMP[ 02 ];  x01i := _TEMP[ 03 ];
    x02r := _TEMP[ 04 ];  x02i := _TEMP[ 05 ];
    x03r := _TEMP[ 06 ];  x03i := _TEMP[ 07 ];
    x04r := _TEMP[ 08 ];  x04i := _TEMP[ 09 ];
    x05r := _TEMP[ 10 ];  x05i := _TEMP[ 11 ];
    {-----------------}   {-----------------}
    x07r := _TEMP[ 14 ];  x07i := _TEMP[ 15 ];
    x08r := _TEMP[ 16 ];  x08i := _TEMP[ 17 ];
    {-----------------}   {-----------------}
    x10r := _TEMP[ 20 ];  x10i := _TEMP[ 21 ];
    x11r := _TEMP[ 22 ];  x11i := _TEMP[ 23 ];
    x12r := _TEMP[ 24 ];  x12i := _TEMP[ 25 ];
    x13r := _TEMP[ 26 ];  x13i := _TEMP[ 27 ];
    x14r := _TEMP[ 28 ];  x14i := _TEMP[ 29 ];
    {-----------------}   {-----------------}

    _TEMP[ 02 ] := x08r;  _TEMP[ 03 ] := x08i;
    _TEMP[ 04 ] := x04r;  _TEMP[ 05 ] := x04i;
    _TEMP[ 06 ] := x12r;  _TEMP[ 07 ] := x12i;
    _TEMP[ 08 ] := x02r;  _TEMP[ 09 ] := x02i;
    _TEMP[ 10 ] := x10r;  _TEMP[ 11 ] := x10i;
    {-----------------}   {-----------------}
    _TEMP[ 14 ] := x14r;  _TEMP[ 15 ] := x14i;
    _TEMP[ 16 ] := x01r;  _TEMP[ 17 ] := x01i;
    {-----------------}   {-----------------}
    _TEMP[ 20 ] := x05r;  _TEMP[ 21 ] := x05i;
    _TEMP[ 22 ] := x13r;  _TEMP[ 23 ] := x13i;
    _TEMP[ 24 ] := x03r;  _TEMP[ 25 ] := x03i;
    _TEMP[ 26 ] := x11r;  _TEMP[ 27 ] := x11i;
    _TEMP[ 28 ] := x07r;  _TEMP[ 29 ] := x07i;
    {-----------------}   {-----------------}
end;

procedure TDiscreteTrans.bitrv216neg;
var
   x01r, x01i,
   x02r, x02i,
   x03r, x03i,
   x04r, x04i,
   x05r, x05i,
   x06r, x06i,
   x07r, x07i,
   x08r, x08i,
   x09r, x09i,
   x10r, x10i,
   x11r, x11i,
   x12r, x12i,
   x13r, x13i,
   x14r, x14i,
   x15r, x15i :Double;
begin
     x01r := _TEMP[ 02 ];  x01i := _TEMP[ 03 ];
     x02r := _TEMP[ 04 ];  x02i := _TEMP[ 05 ];
     x03r := _TEMP[ 06 ];  x03i := _TEMP[ 07 ];
     x04r := _TEMP[ 08 ];  x04i := _TEMP[ 09 ];
     x05r := _TEMP[ 10 ];  x05i := _TEMP[ 11 ];
     x06r := _TEMP[ 12 ];  x06i := _TEMP[ 13 ];
     x07r := _TEMP[ 14 ];  x07i := _TEMP[ 15 ];
     x08r := _TEMP[ 16 ];  x08i := _TEMP[ 17 ];
     x09r := _TEMP[ 18 ];  x09i := _TEMP[ 19 ];
     x10r := _TEMP[ 20 ];  x10i := _TEMP[ 21 ];
     x11r := _TEMP[ 22 ];  x11i := _TEMP[ 23 ];
     x12r := _TEMP[ 24 ];  x12i := _TEMP[ 25 ];
     x13r := _TEMP[ 26 ];  x13i := _TEMP[ 27 ];
     x14r := _TEMP[ 28 ];  x14i := _TEMP[ 29 ];
     x15r := _TEMP[ 30 ];  x15i := _TEMP[ 31 ];

     _TEMP[ 02 ] := x15r;  _TEMP[ 03 ] := x15i;
     _TEMP[ 04 ] := x07r;  _TEMP[ 05 ] := x07i;
     _TEMP[ 06 ] := x11r;  _TEMP[ 07 ] := x11i;
     _TEMP[ 08 ] := x03r;  _TEMP[ 09 ] := x03i;
     _TEMP[ 10 ] := x13r;  _TEMP[ 11 ] := x13i;
     _TEMP[ 12 ] := x05r;  _TEMP[ 13 ] := x05i;
     _TEMP[ 14 ] := x09r;  _TEMP[ 15 ] := x09i;
     _TEMP[ 16 ] := x01r;  _TEMP[ 17 ] := x01i;
     _TEMP[ 18 ] := x14r;  _TEMP[ 19 ] := x14i;
     _TEMP[ 20 ] := x06r;  _TEMP[ 21 ] := x06i;
     _TEMP[ 22 ] := x10r;  _TEMP[ 23 ] := x10i;
     _TEMP[ 24 ] := x02r;  _TEMP[ 25 ] := x02i;
     _TEMP[ 26 ] := x12r;  _TEMP[ 27 ] := x12i;
     _TEMP[ 28 ] := x04r;  _TEMP[ 29 ] := x04i;
     _TEMP[ 30 ] := x08r;  _TEMP[ 31 ] := x08i;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.bitrv208;
var
   x1r, x1i,
   {-}  {-}
   x3r, x3i,
   x4r, x4i,
   {-}  {-}
   x6r, x6i :Double;
begin
     x1r := _TEMP[ 02 ];  x1i := _TEMP[ 03 ];
     {----------------}   {----------------}
     x3r := _TEMP[ 06 ];  x3i := _TEMP[ 07 ];
     x4r := _TEMP[ 08 ];  x4i := _TEMP[ 09 ];
     {----------------}   {----------------}
     x6r := _TEMP[ 12 ];  x6i := _TEMP[ 13 ];

     _TEMP[ 02 ] := x4r;  _TEMP[ 03 ] := x4i;
     {----------------}   {----------------}
     _TEMP[ 06 ] := x6r;  _TEMP[ 07 ] := x6i;
     _TEMP[ 08 ] := x1r;  _TEMP[ 09 ] := x1i;
     {----------------}   {----------------}
     _TEMP[ 12 ] := x3r;  _TEMP[ 13 ] := x3i;
end;

procedure TDiscreteTrans.bitrv208neg;
var
   x1r, x1i,
   x2r, x2i,
   x3r, x3i,
   x4r, x4i,
   x5r, x5i,
   x6r, x6i,
   x7r, x7i :Double;
begin
     x1r := _TEMP[ 02 ];  x1i := _TEMP[ 03 ];
     x2r := _TEMP[ 04 ];  x2i := _TEMP[ 05 ];
     x3r := _TEMP[ 06 ];  x3i := _TEMP[ 07 ];
     x4r := _TEMP[ 08 ];  x4i := _TEMP[ 09 ];
     x5r := _TEMP[ 10 ];  x5i := _TEMP[ 11 ];
     x6r := _TEMP[ 12 ];  x6i := _TEMP[ 13 ];
     x7r := _TEMP[ 14 ];  x7i := _TEMP[ 15 ];

     _TEMP[ 02 ] := x7r;  _TEMP[ 03 ] := x7i;
     _TEMP[ 04 ] := x3r;  _TEMP[ 05 ] := x3i;
     _TEMP[ 06 ] := x5r;  _TEMP[ 07 ] := x5i;
     _TEMP[ 08 ] := x1r;  _TEMP[ 09 ] := x1i;
     _TEMP[ 10 ] := x6r;  _TEMP[ 11 ] := x6i;
     _TEMP[ 12 ] := x2r;  _TEMP[ 13 ] := x2i;
     _TEMP[ 14 ] := x4r;  _TEMP[ 15 ] := x4i;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftf1st( var w:array of Double );
var
   j, j0, j1, j2, j3, k, m, mh :Integer;
   wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i,
   wd1r, wd1i, wd3r, wd3i,
   x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
   y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i :Double;
begin
     mh := _TempN shr 3;
     m := 2 * mh;

     j1 :=      m;
     j2 := j1 + m;
     j3 := j2 + m;

     x0r := _TEMP[ 0  ] + _TEMP[ j2 ]; x0i := _TEMP[ 1      ] + _TEMP[ j2 + 1 ];
     x1r := _TEMP[ 0  ] - _TEMP[ j2 ]; x1i := _TEMP[ 1      ] - _TEMP[ j2 + 1 ];
     x2r := _TEMP[ j1 ] + _TEMP[ j3 ]; x2i := _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
     x3r := _TEMP[ j1 ] - _TEMP[ j3 ]; x3i := _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
     _TEMP[ 0 ] := x0r + x2r;
     _TEMP[ 1 ] := x0i + x2i;
     _TEMP[ j1 ] := x0r - x2r;
     _TEMP[ j1 + 1 ] := x0i - x2i;
     _TEMP[ j2 ] := x1r - x3i;
     _TEMP[ j2 + 1 ] := x1i + x3r;
     _TEMP[ j3 ] := x1r + x3i;
     _TEMP[ j3 + 1 ] := x1i - x3r;

     wn4r := w[ 1 ];
     csc1 := w[ 2 ];
     csc3 := w[ 3 ];
     wd1r := 1;
     wd1i := 0;
     wd3r := 1;
     wd3i := 0;
     k := 0;
     j := 2;
     while j < mh - 2 do
     begin
          Inc( k, 4 );

          wk1r := csc1 * ( wd1r + w[ k     ] ); wk1i := csc1 * ( wd1i + w[ k + 1 ] );
          wk3r := csc3 * ( wd3r + w[ k + 2 ] ); wk3i := csc3 * ( wd3i + w[ k + 3 ] );

          wd1r := w[ k     ]; wd1i := w[ k + 1 ];
          wd3r := w[ k + 2 ]; wd3i := w[ k + 3 ];

          j1 := j  + m;
          j2 := j1 + m;
          j3 := j2 + m;

          x0r := _TEMP[ j      ] + _TEMP[ j2     ]; x0i := _TEMP[ j  + 1 ] + _TEMP[ j2 + 1 ];
          x1r := _TEMP[ j      ] - _TEMP[ j2     ]; x1i := _TEMP[ j  + 1 ] - _TEMP[ j2 + 1 ];
          y0r := _TEMP[ j  + 2 ] + _TEMP[ j2 + 2 ]; y0i := _TEMP[ j  + 3 ] + _TEMP[ j2 + 3 ];
          y1r := _TEMP[ j  + 2 ] - _TEMP[ j2 + 2 ]; y1i := _TEMP[ j  + 3 ] - _TEMP[ j2 + 3 ];
          x2r := _TEMP[ j1     ] + _TEMP[ j3     ]; x2i := _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
          x3r := _TEMP[ j1     ] - _TEMP[ j3     ]; x3i := _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
          y2r := _TEMP[ j1 + 2 ] + _TEMP[ j3 + 2 ]; y2i := _TEMP[ j1 + 3 ] + _TEMP[ j3 + 3 ];
          y3r := _TEMP[ j1 + 2 ] - _TEMP[ j3 + 2 ]; y3i := _TEMP[ j1 + 3 ] - _TEMP[ j3 + 3 ];
          _TEMP[ j ] := x0r + x2r;
          _TEMP[ j + 1 ] := x0i + x2i;
          _TEMP[ j + 2 ] := y0r + y2r;
          _TEMP[ j + 3 ] := y0i + y2i;
          _TEMP[ j1 ] := x0r - x2r;
          _TEMP[ j1 + 1 ] := x0i - x2i;
          _TEMP[ j1 + 2 ] := y0r - y2r;
          _TEMP[ j1 + 3 ] := y0i - y2i;

          x0r := x1r - x3i; x0i := x1i + x3r;
          _TEMP[ j2     ] := wk1r * x0r - wk1i * x0i;
          _TEMP[ j2 + 1 ] := wk1r * x0i + wk1i * x0r;

          x0r := y1r - y3i; x0i := y1i + y3r;
          _TEMP[ j2 + 2 ] := wd1r * x0r - wd1i * x0i;
          _TEMP[ j2 + 3 ] := wd1r * x0i + wd1i * x0r;

          x0r := x1r + x3i; x0i := x1i - x3r;
          _TEMP[ j3 ] := wk3r * x0r + wk3i * x0i;
          _TEMP[ j3 + 1 ] := wk3r * x0i - wk3i * x0r;

          x0r := y1r + y3i; x0i := y1i - y3r;
          _TEMP[ j3 + 2 ] := wd3r * x0r + wd3i * x0i;
          _TEMP[ j3 + 3 ] := wd3r * x0i - wd3i * x0r;

          j0 :=  m - j;
          j1 := j0 + m;
          j2 := j1 + m;
          j3 := j2 + m;

          x0r := _TEMP[ j0     ] + _TEMP[ j2     ]; x0i := _TEMP[ j0 + 1 ] + _TEMP[ j2 + 1 ];
          x1r := _TEMP[ j0     ] - _TEMP[ j2     ]; x1i := _TEMP[ j0 + 1 ] - _TEMP[ j2 + 1 ];
          y0r := _TEMP[ j0 - 2 ] + _TEMP[ j2 - 2 ]; y0i := _TEMP[ j0 - 1 ] + _TEMP[ j2 - 1 ];
          y1r := _TEMP[ j0 - 2 ] - _TEMP[ j2 - 2 ]; y1i := _TEMP[ j0 - 1 ] - _TEMP[ j2 - 1 ];
          x2r := _TEMP[ j1     ] + _TEMP[ j3     ]; x2i := _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
          x3r := _TEMP[ j1     ] - _TEMP[ j3     ]; x3i := _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
          y2r := _TEMP[ j1 - 2 ] + _TEMP[ j3 - 2 ]; y2i := _TEMP[ j1 - 1 ] + _TEMP[ j3 - 1 ];
          y3r := _TEMP[ j1 - 2 ] - _TEMP[ j3 - 2 ]; y3i := _TEMP[ j1 - 1 ] - _TEMP[ j3 - 1 ];
          _TEMP[ j0     ] := x0r + x2r;
          _TEMP[ j0 + 1 ] := x0i + x2i;
          _TEMP[ j0 - 2 ] := y0r + y2r;
          _TEMP[ j0 - 1 ] := y0i + y2i;
          _TEMP[ j1     ] := x0r - x2r;
          _TEMP[ j1 + 1 ] := x0i - x2i;
          _TEMP[ j1 - 2 ] := y0r - y2r;
          _TEMP[ j1 - 1 ] := y0i - y2i;

          x0r := x1r - x3i; x0i := x1i + x3r;
          _TEMP[ j2     ] := wk1i * x0r - wk1r * x0i;
          _TEMP[ j2 + 1 ] := wk1i * x0i + wk1r * x0r;

          x0r := y1r - y3i; x0i := y1i + y3r;
          _TEMP[ j2 - 2 ] := wd1i * x0r - wd1r * x0i;
          _TEMP[ j2 - 1 ] := wd1i * x0i + wd1r * x0r;

          x0r := x1r + x3i; x0i := x1i - x3r;
          _TEMP[ j3     ] := wk3i * x0r + wk3r * x0i;
          _TEMP[ j3 + 1 ] := wk3i * x0i - wk3r * x0r;

          x0r := y1r + y3i; x0i := y1i - y3r;
          _TEMP[ j3 - 2 ] := wd3i * x0r + wd3r * x0i;
          _TEMP[ j3 - 1 ] := wd3i * x0i - wd3r * x0r;

          Inc( j, 4 )
     end;

     wk1r := csc1 * ( wd1r + wn4r );
     wk1i := csc1 * ( wd1i + wn4r );
     wk3r := csc3 * ( wd3r - wn4r );
     wk3i := csc3 * ( wd3i - wn4r );

     j0 := mh;
     j1 := j0 + m;
     j2 := j1 + m;
     j3 := j2 + m;

     x0r := _TEMP[ j0 - 2 ] + _TEMP[ j2 - 2 ]; x0i := _TEMP[ j0 - 1 ] + _TEMP[ j2 - 1 ];
     x1r := _TEMP[ j0 - 2 ] - _TEMP[ j2 - 2 ]; x1i := _TEMP[ j0 - 1 ] - _TEMP[ j2 - 1 ];
     x2r := _TEMP[ j1 - 2 ] + _TEMP[ j3 - 2 ]; x2i := _TEMP[ j1 - 1 ] + _TEMP[ j3 - 1 ];
     x3r := _TEMP[ j1 - 2 ] - _TEMP[ j3 - 2 ]; x3i := _TEMP[ j1 - 1 ] - _TEMP[ j3 - 1 ];
     _TEMP[ j0 - 2 ] := x0r + x2r;
     _TEMP[ j0 - 1 ] := x0i + x2i;
     _TEMP[ j1 - 2 ] := x0r - x2r;
     _TEMP[ j1 - 1 ] := x0i - x2i;

     x0r := x1r - x3i; x0i := x1i + x3r;
     _TEMP[ j2 - 2 ] := wk1r * x0r - wk1i * x0i;
     _TEMP[ j2 - 1 ] := wk1r * x0i + wk1i * x0r;

     x0r := x1r + x3i; x0i := x1i - x3r;
     _TEMP[ j3 - 2 ] := wk3r * x0r + wk3i * x0i;
     _TEMP[ j3 - 1 ] := wk3r * x0i - wk3i * x0r;

     x0r := _TEMP[ j0     ] + _TEMP[ j2     ]; x0i := _TEMP[ j0 + 1 ] + _TEMP[ j2 + 1 ];
     x1r := _TEMP[ j0     ] - _TEMP[ j2     ]; x1i := _TEMP[ j0 + 1 ] - _TEMP[ j2 + 1 ];
     x2r := _TEMP[ j1     ] + _TEMP[ j3     ]; x2i := _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
     x3r := _TEMP[ j1     ] - _TEMP[ j3     ]; x3i := _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
     _TEMP[ j0     ] := x0r + x2r;
     _TEMP[ j0 + 1 ] := x0i + x2i;
     _TEMP[ j1     ] := x0r - x2r;
     _TEMP[ j1 + 1 ] := x0i - x2i;

     x0r := x1r - x3i; x0i := x1i + x3r;
     _TEMP[ j2     ] := wn4r * ( x0r - x0i );
     _TEMP[ j2 + 1 ] := wn4r * ( x0i + x0r );

     x0r := x1r + x3i; x0i := x1i - x3r;
     _TEMP[ j3     ] := -wn4r * ( x0r + x0i );
     _TEMP[ j3 + 1 ] := -wn4r * ( x0i - x0r );

     x0r := _TEMP[ j0 + 2 ] + _TEMP[ j2 + 2 ]; x0i := _TEMP[ j0 + 3 ] + _TEMP[ j2 + 3 ];
     x1r := _TEMP[ j0 + 2 ] - _TEMP[ j2 + 2 ]; x1i := _TEMP[ j0 + 3 ] - _TEMP[ j2 + 3 ];
     x2r := _TEMP[ j1 + 2 ] + _TEMP[ j3 + 2 ]; x2i := _TEMP[ j1 + 3 ] + _TEMP[ j3 + 3 ];
     x3r := _TEMP[ j1 + 2 ] - _TEMP[ j3 + 2 ]; x3i := _TEMP[ j1 + 3 ] - _TEMP[ j3 + 3 ];
     _TEMP[ j0 + 2 ] := x0r + x2r;
     _TEMP[ j0 + 3 ] := x0i + x2i;
     _TEMP[ j1 + 2 ] := x0r - x2r;
     _TEMP[ j1 + 3 ] := x0i - x2i;

     x0r := x1r - x3i; x0i := x1i + x3r;
     _TEMP[ j2 + 2 ] := wk1i * x0r - wk1r * x0i;
     _TEMP[ j2 + 3 ] := wk1i * x0i + wk1r * x0r;

     x0r := x1r + x3i; x0i := x1i - x3r;
     _TEMP[ j3 + 2 ] := wk3i * x0r + wk3r * x0i;
     _TEMP[ j3 + 3 ] := wk3i * x0i - wk3r * x0r;
end;

procedure TDiscreteTrans.cftb1st( var w:array of Double );
var
   j, j0, j1, j2, j3, k, m, mh :Integer;
   wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i,
   wd1r, wd1i, wd3r, wd3i,
   x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
   y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i :Double;
begin
     mh := _TempN shr 3;
     m := 2 * mh;

     j1 := m;
     j2 := j1 + m;
     j3 := j2 + m;

     x0r :=  _TEMP[      0 ] + _TEMP[ j2     ]; x0i := -_TEMP[      1 ] - _TEMP[ j2 + 1 ];
     x1r :=  _TEMP[      0 ] - _TEMP[ j2     ]; x1i := -_TEMP[      1 ] + _TEMP[ j2 + 1 ];
     x2r :=  _TEMP[ j1     ] + _TEMP[ j3     ]; x2i :=  _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
     x3r :=  _TEMP[ j1     ] - _TEMP[ j3     ]; x3i :=  _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
     _TEMP[      0 ] := x0r + x2r;
     _TEMP[      1 ] := x0i - x2i;
     _TEMP[ j1     ] := x0r - x2r;
     _TEMP[ j1 + 1 ] := x0i + x2i;
     _TEMP[ j2     ] := x1r + x3i;
     _TEMP[ j2 + 1 ] := x1i + x3r;
     _TEMP[ j3     ] := x1r - x3i;
     _TEMP[ j3 + 1 ] := x1i - x3r;

     wn4r := w[ 1 ];
     csc1 := w[ 2 ];
     csc3 := w[ 3 ];
     wd1r := 1;
     wd1i := 0;
     wd3r := 1;
     wd3i := 0;
     k := 0;
     j := 2;
     while j < mh - 2 do
     begin
          Inc( k, 4 );

          wk1r := csc1 * ( wd1r + w[ k     ] ); wk1i := csc1 * ( wd1i + w[ k + 1 ] );
          wk3r := csc3 * ( wd3r + w[ k + 2 ] ); wk3i := csc3 * ( wd3i + w[ k + 3 ] );
          wd1r := w[ k     ]; wd1i := w[ k + 1 ];
          wd3r := w[ k + 2 ]; wd3i := w[ k + 3 ];

          j1 := j  + m;
          j2 := j1 + m;
          j3 := j2 + m;

          x0r :=  _TEMP[ j      ] + _TEMP[ j2     ];
          x0i := -_TEMP[ j  + 1 ] - _TEMP[ j2 + 1 ];
          x1r :=  _TEMP[ j      ] - _TEMP[ j2     ];
          x1i := -_TEMP[ j  + 1 ] + _TEMP[ j2 + 1 ];
          y0r :=  _TEMP[ j  + 2 ] + _TEMP[ j2 + 2 ];
          y0i := -_TEMP[ j  + 3 ] - _TEMP[ j2 + 3 ];
          y1r :=  _TEMP[ j  + 2 ] - _TEMP[ j2 + 2 ];
          y1i := -_TEMP[ j  + 3 ] + _TEMP[ j2 + 3 ];
          x2r :=  _TEMP[ j1     ] + _TEMP[ j3     ];
          x2i :=  _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
          x3r :=  _TEMP[ j1     ] - _TEMP[ j3     ];
          x3i :=  _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
          y2r :=  _TEMP[ j1 + 2 ] + _TEMP[ j3 + 2 ];
          y2i :=  _TEMP[ j1 + 3 ] + _TEMP[ j3 + 3 ];
          y3r :=  _TEMP[ j1 + 2 ] - _TEMP[ j3 + 2 ];
          y3i :=  _TEMP[ j1 + 3 ] - _TEMP[ j3 + 3 ];
          _TEMP[ j      ] := x0r + x2r;
          _TEMP[ j  + 1 ] := x0i - x2i;
          _TEMP[ j  + 2 ] := y0r + y2r;
          _TEMP[ j  + 3 ] := y0i - y2i;
          _TEMP[ j1     ] := x0r - x2r;
          _TEMP[ j1 + 1 ] := x0i + x2i;
          _TEMP[ j1 + 2 ] := y0r - y2r;
          _TEMP[ j1 + 3 ] := y0i + y2i;

          x0r := x1r + x3i; x0i := x1i + x3r;
          _TEMP[ j2     ] := wk1r * x0r - wk1i * x0i;
          _TEMP[ j2 + 1 ] := wk1r * x0i + wk1i * x0r;

          x0r := y1r + y3i; x0i := y1i + y3r;
          _TEMP[ j2 + 2 ] := wd1r * x0r - wd1i * x0i;
          _TEMP[ j2 + 3 ] := wd1r * x0i + wd1i * x0r;

          x0r := x1r - x3i; x0i := x1i - x3r;
          _TEMP[ j3     ] := wk3r * x0r + wk3i * x0i;
          _TEMP[ j3 + 1 ] := wk3r * x0i - wk3i * x0r;

          x0r := y1r - y3i; x0i := y1i - y3r;
          _TEMP[ j3 + 2 ] := wd3r * x0r + wd3i * x0i;
          _TEMP[ j3 + 3 ] := wd3r * x0i - wd3i * x0r;

          j0 :=  m - j;
          j1 := j0 + m;
          j2 := j1 + m;
          j3 := j2 + m;

          x0r :=  _TEMP[ j0     ] + _TEMP[ j2     ]; x0i := -_TEMP[ j0 + 1 ] - _TEMP[ j2 + 1 ];
          x1r :=  _TEMP[ j0     ] - _TEMP[ j2     ]; x1i := -_TEMP[ j0 + 1 ] + _TEMP[ j2 + 1 ];
          y0r :=  _TEMP[ j0 - 2 ] + _TEMP[ j2 - 2 ]; y0i := -_TEMP[ j0 - 1 ] - _TEMP[ j2 - 1 ];
          y1r :=  _TEMP[ j0 - 2 ] - _TEMP[ j2 - 2 ]; y1i := -_TEMP[ j0 - 1 ] + _TEMP[ j2 - 1 ];
          x2r :=  _TEMP[ j1     ] + _TEMP[ j3     ]; x2i :=  _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
          x3r :=  _TEMP[ j1     ] - _TEMP[ j3     ]; x3i :=  _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
          y2r :=  _TEMP[ j1 - 2 ] + _TEMP[ j3 - 2 ]; y2i :=  _TEMP[ j1 - 1 ] + _TEMP[ j3 - 1 ];
          y3r :=  _TEMP[ j1 - 2 ] - _TEMP[ j3 - 2 ]; y3i :=  _TEMP[ j1 - 1 ] - _TEMP[ j3 - 1 ];
          _TEMP[ j0     ] := x0r + x2r;
          _TEMP[ j0 + 1 ] := x0i - x2i;
          _TEMP[ j0 - 2 ] := y0r + y2r;
          _TEMP[ j0 - 1 ] := y0i - y2i;
          _TEMP[ j1     ] := x0r - x2r;
          _TEMP[ j1 + 1 ] := x0i + x2i;
          _TEMP[ j1 - 2 ] := y0r - y2r;
          _TEMP[ j1 - 1 ] := y0i + y2i;

          x0r := x1r + x3i; x0i := x1i + x3r;
          _TEMP[ j2 ] := wk1i * x0r - wk1r * x0i;
          _TEMP[ j2 + 1 ] := wk1i * x0i + wk1r * x0r;

          x0r := y1r + y3i; x0i := y1i + y3r;
          _TEMP[ j2 - 2 ] := wd1i * x0r - wd1r * x0i;
          _TEMP[ j2 - 1 ] := wd1i * x0i + wd1r * x0r;

          x0r := x1r - x3i; x0i := x1i - x3r;
          _TEMP[ j3 ] := wk3i * x0r + wk3r * x0i;
          _TEMP[ j3 + 1 ] := wk3i * x0i - wk3r * x0r;

          x0r := y1r - y3i; x0i := y1i - y3r;
          _TEMP[ j3 - 2 ] := wd3i * x0r + wd3r * x0i;
          _TEMP[ j3 - 1 ] := wd3i * x0i - wd3r * x0r;

          Inc( j, 4 )
     end;

     wk1r := csc1 * ( wd1r + wn4r );
     wk1i := csc1 * ( wd1i + wn4r );
     wk3r := csc3 * ( wd3r - wn4r );
     wk3i := csc3 * ( wd3i - wn4r );

     j0 := mh;
     j1 := j0 + m;
     j2 := j1 + m;
     j3 := j2 + m;

     x0r :=  _TEMP[ j0 - 2 ] + _TEMP[ j2 - 2 ];  x0i := -_TEMP[ j0 - 1 ] - _TEMP[ j2 - 1 ];
     x1r :=  _TEMP[ j0 - 2 ] - _TEMP[ j2 - 2 ];  x1i := -_TEMP[ j0 - 1 ] + _TEMP[ j2 - 1 ];
     x2r :=  _TEMP[ j1 - 2 ] + _TEMP[ j3 - 2 ];  x2i :=  _TEMP[ j1 - 1 ] + _TEMP[ j3 - 1 ];
     x3r :=  _TEMP[ j1 - 2 ] - _TEMP[ j3 - 2 ];  x3i :=  _TEMP[ j1 - 1 ] - _TEMP[ j3 - 1 ];
     _TEMP[ j0 - 2 ] := x0r + x2r;
     _TEMP[ j0 - 1 ] := x0i - x2i;
     _TEMP[ j1 - 2 ] := x0r - x2r;
     _TEMP[ j1 - 1 ] := x0i + x2i;

     x0r := x1r + x3i; x0i := x1i + x3r;
     _TEMP[ j2 - 2 ] := wk1r * x0r - wk1i * x0i;
     _TEMP[ j2 - 1 ] := wk1r * x0i + wk1i * x0r;

     x0r := x1r - x3i; x0i := x1i - x3r;
     _TEMP[ j3 - 2 ] := wk3r * x0r + wk3i * x0i;
     _TEMP[ j3 - 1 ] := wk3r * x0i - wk3i * x0r;

     x0r :=  _TEMP[ j0     ] + _TEMP[ j2     ];  x0i := -_TEMP[ j0 + 1 ] - _TEMP[ j2 + 1 ];
     x1r :=  _TEMP[ j0     ] - _TEMP[ j2     ];  x1i := -_TEMP[ j0 + 1 ] + _TEMP[ j2 + 1 ];
     x2r :=  _TEMP[ j1     ] + _TEMP[ j3     ];  x2i :=  _TEMP[ j1 + 1 ] + _TEMP[ j3 + 1 ];
     x3r :=  _TEMP[ j1     ] - _TEMP[ j3     ];  x3i :=  _TEMP[ j1 + 1 ] - _TEMP[ j3 + 1 ];
     _TEMP[ j0     ] := x0r + x2r;
     _TEMP[ j0 + 1 ] := x0i - x2i;
     _TEMP[ j1     ] := x0r - x2r;
     _TEMP[ j1 + 1 ] := x0i + x2i;

     x0r := x1r + x3i; x0i := x1i + x3r;
     _TEMP[ j2 ] := wn4r * ( x0r - x0i );
     _TEMP[ j2 + 1 ] := wn4r * ( x0i + x0r );

     x0r := x1r - x3i; x0i := x1i - x3r;
     _TEMP[ j3 ] := -wn4r * ( x0r + x0i );
     _TEMP[ j3 + 1 ] := -wn4r * ( x0i - x0r );

     x0r :=  _TEMP[ j0 + 2 ] + _TEMP[ j2 + 2 ];  x0i := -_TEMP[ j0 + 3 ] - _TEMP[ j2 + 3 ];
     x1r :=  _TEMP[ j0 + 2 ] - _TEMP[ j2 + 2 ];  x1i := -_TEMP[ j0 + 3 ] + _TEMP[ j2 + 3 ];
     x2r :=  _TEMP[ j1 + 2 ] + _TEMP[ j3 + 2 ];  x2i :=  _TEMP[ j1 + 3 ] + _TEMP[ j3 + 3 ];
     x3r :=  _TEMP[ j1 + 2 ] - _TEMP[ j3 + 2 ];  x3i :=  _TEMP[ j1 + 3 ] - _TEMP[ j3 + 3 ];
     _TEMP[ j0 + 2 ] := x0r + x2r;
     _TEMP[ j0 + 3 ] := x0i - x2i;
     _TEMP[ j1 + 2 ] := x0r - x2r;
     _TEMP[ j1 + 3 ] := x0i + x2i;

     x0r := x1r + x3i; x0i := x1i + x3r;
     _TEMP[ j2 + 2 ] := wk1i * x0r - wk1r * x0i;
     _TEMP[ j2 + 3 ] := wk1i * x0i + wk1r * x0r;

     x0r := x1r - x3i; x0i := x1i - x3r;
     _TEMP[ j3 + 2 ] := wk3i * x0r + wk3r * x0i;
     _TEMP[ j3 + 3 ] := wk3i * x0i - wk3r * x0r;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftmdl1( const n:Integer; var a:array of Double; var w:array of Double );
var
   j, j0, j1, j2, j3, k, m, mh :Integer;
   wn4r, wk1r, wk1i, wk3r, wk3i,
   x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i :Double;
begin
    mh := n shr 3;
    m := 2 * mh;

    j1 :=      m;
    j2 := j1 + m;
    j3 := j2 + m;

    x0r := a[      0 ] + a[ j2     ];
    x0i := a[      1 ] + a[ j2 + 1 ];
    x1r := a[      0 ] - a[ j2     ];
    x1i := a[      1 ] - a[ j2 + 1 ];
    x2r := a[ j1     ] + a[ j3     ];
    x2i := a[ j1 + 1 ] + a[ j3 + 1 ];
    x3r := a[ j1     ] - a[ j3     ];
    x3i := a[ j1 + 1 ] - a[ j3 + 1 ];
    a[      0 ] := x0r + x2r;
    a[      1 ] := x0i + x2i;
    a[ j1     ] := x0r - x2r;
    a[ j1 + 1 ] := x0i - x2i;
    a[ j2     ] := x1r - x3i;
    a[ j2 + 1 ] := x1i + x3r;
    a[ j3     ] := x1r + x3i;
    a[ j3 + 1 ] := x1i - x3r;

    wn4r := w[ 1 ];
    k := 0;
    j := 2;
    while j < mh do
    begin
        Inc( k, 4 );

        wk1r := w[ k     ];
        wk1i := w[ k + 1 ];
        wk3r := w[ k + 2 ];
        wk3i := w[ k + 3 ];

        j1 := j  + m;
        j2 := j1 + m;
        j3 := j2 + m;

        x0r := a[ j      ] + a[ j2     ];
        x0i := a[ j  + 1 ] + a[ j2 + 1 ];
        x1r := a[ j      ] - a[ j2     ];
        x1i := a[ j  + 1 ] - a[ j2 + 1 ];
        x2r := a[ j1     ] + a[ j3     ];
        x2i := a[ j1 + 1 ] + a[ j3 + 1 ];
        x3r := a[ j1     ] - a[ j3     ];
        x3i := a[ j1 + 1 ] - a[ j3 + 1 ];
        a[ j      ] := x0r + x2r;
        a[ j  + 1 ] := x0i + x2i;
        a[ j1     ] := x0r - x2r;
        a[ j1 + 1 ] := x0i - x2i;

        x0r := x1r - x3i;
        x0i := x1i + x3r;
        a[ j2     ] := wk1r * x0r - wk1i * x0i;
        a[ j2 + 1 ] := wk1r * x0i + wk1i * x0r;

        x0r := x1r + x3i;
        x0i := x1i - x3r;
        a[ j3     ] := wk3r * x0r + wk3i * x0i;
        a[ j3 + 1 ] := wk3r * x0i - wk3i * x0r;

        j0 := m  - j;
        j1 := j0 + m;
        j2 := j1 + m;
        j3 := j2 + m;

        x0r := a[ j0     ] + a[ j2     ];
        x0i := a[ j0 + 1 ] + a[ j2 + 1 ];
        x1r := a[ j0     ] - a[ j2     ];
        x1i := a[ j0 + 1 ] - a[ j2 + 1 ];
        x2r := a[ j1     ] + a[ j3     ];
        x2i := a[ j1 + 1 ] + a[ j3 + 1 ];
        x3r := a[ j1     ] - a[ j3     ];
        x3i := a[ j1 + 1 ] - a[ j3 + 1 ];
        a[ j0     ] := x0r + x2r;
        a[ j0 + 1 ] := x0i + x2i;
        a[ j1     ] := x0r - x2r;
        a[ j1 + 1 ] := x0i - x2i;

        x0r := x1r - x3i;
        x0i := x1i + x3r;
        a[ j2     ] := wk1i * x0r - wk1r * x0i;
        a[ j2 + 1 ] := wk1i * x0i + wk1r * x0r;

        x0r := x1r + x3i;
        x0i := x1i - x3r;
        a[ j3     ] := wk3i * x0r + wk3r * x0i;
        a[ j3 + 1 ] := wk3i * x0i - wk3r * x0r;

        Inc( j, 2 )
    end;

    j0 := mh;
    j1 := j0 + m;
    j2 := j1 + m;
    j3 := j2 + m;

    x0r := a[ j0     ] + a[ j2     ];
    x0i := a[ j0 + 1 ] + a[ j2 + 1 ];
    x1r := a[ j0     ] - a[ j2     ];
    x1i := a[ j0 + 1 ] - a[ j2 + 1 ];
    x2r := a[ j1     ] + a[ j3     ];
    x2i := a[ j1 + 1 ] + a[ j3 + 1 ];
    x3r := a[ j1     ] - a[ j3     ];
    x3i := a[ j1 + 1 ] - a[ j3 + 1 ];
    a[ j0     ] := x0r + x2r;
    a[ j0 + 1 ] := x0i + x2i;
    a[ j1     ] := x0r - x2r;
    a[ j1 + 1 ] := x0i - x2i;

    x0r := x1r - x3i;
    x0i := x1i + x3r;
    a[ j2     ] := +wn4r * ( x0r - x0i );
    a[ j2 + 1 ] := +wn4r * ( x0i + x0r );

    x0r := x1r + x3i;
    x0i := x1i - x3r;
    a[ j3     ] := -wn4r * ( x0r + x0i );
    a[ j3 + 1 ] := -wn4r * ( x0i - x0r );
end;

procedure TDiscreteTrans.cftmdl2( const n:Integer; var a:array of Double; var w:array of Double );
var
   j, j0, j1, j2, j3, k, kr, m, mh :Integer;
   wn4r, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i,
   x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i, y0r, y0i, y2r, y2i :Double;
 begin
    mh := n shr 3;
    m := 2 * mh;

    wn4r := w[ 1 ];

    j1 := m;
    j2 := j1 + m;
    j3 := j2 + m;

    x0r := a[      0 ] - a[ j2 + 1 ];
    x0i := a[      1 ] + a[ j2     ];
    x1r := a[      0 ] + a[ j2 + 1 ];
    x1i := a[      1 ] - a[ j2     ];
    x2r := a[ j1     ] - a[ j3 + 1 ];
    x2i := a[ j1 + 1 ] + a[ j3     ];
    x3r := a[ j1     ] + a[ j3 + 1 ];
    x3i := a[ j1 + 1 ] - a[ j3     ];
    y0r := wn4r * ( x2r - x2i );
    y0i := wn4r * ( x2i + x2r );
    a[      0 ] := x0r + y0r;
    a[      1 ] := x0i + y0i;
    a[ j1     ] := x0r - y0r;
    a[ j1 + 1 ] := x0i - y0i;

    y0r := wn4r * ( x3r - x3i );
    y0i := wn4r * ( x3i + x3r );
    a[ j2     ] := x1r - y0i;
    a[ j2 + 1 ] := x1i + y0r;
    a[ j3     ] := x1r + y0i;
    a[ j3 + 1 ] := x1i - y0r;

    k := 0;
    kr := 2 * m;
    j := 2;
    while j < mh do
    begin
        Inc( k, 4 );

        wk1r := w[ k     ];
        wk1i := w[ k + 1 ];
        wk3r := w[ k + 2 ];
        wk3i := w[ k + 3 ];

        Dec( kr, 4 );
        wd1i := w[ kr     ];
        wd1r := w[ kr + 1 ];
        wd3i := w[ kr + 2 ];
        wd3r := w[ kr + 3 ];

        j1 := j  + m;
        j2 := j1 + m;
        j3 := j2 + m;

        x0r := a[ j      ] - a[ j2 + 1 ];
        x0i := a[ j  + 1 ] + a[ j2     ];
        x1r := a[ j      ] + a[ j2 + 1 ];
        x1i := a[ j  + 1 ] - a[ j2     ];
        x2r := a[ j1     ] - a[ j3 + 1 ];
        x2i := a[ j1 + 1 ] + a[ j3     ];
        x3r := a[ j1     ] + a[ j3 + 1 ];
        x3i := a[ j1 + 1 ] - a[ j3     ];

        y0r := wk1r * x0r - wk1i * x0i;
        y0i := wk1r * x0i + wk1i * x0r;
        y2r := wd1r * x2r - wd1i * x2i;
        y2i := wd1r * x2i + wd1i * x2r;
        a[ j      ] := y0r + y2r;
        a[ j  + 1 ] := y0i + y2i;
        a[ j1     ] := y0r - y2r;
        a[ j1 + 1 ] := y0i - y2i;

        y0r := wk3r * x1r + wk3i * x1i;
        y0i := wk3r * x1i - wk3i * x1r;
        y2r := wd3r * x3r + wd3i * x3i;
        y2i := wd3r * x3i - wd3i * x3r;
        a[ j2     ] := y0r + y2r;
        a[ j2 + 1 ] := y0i + y2i;
        a[ j3     ] := y0r - y2r;
        a[ j3 + 1 ] := y0i - y2i;

        j0 := m  - j;
        j1 := j0 + m;
        j2 := j1 + m;
        j3 := j2 + m;

        x0r := a[ j0     ] - a[ j2 + 1 ];
        x0i := a[ j0 + 1 ] + a[ j2     ];
        x1r := a[ j0     ] + a[ j2 + 1 ];
        x1i := a[ j0 + 1 ] - a[ j2     ];
        x2r := a[ j1     ] - a[ j3 + 1 ];
        x2i := a[ j1 + 1 ] + a[ j3     ];
        x3r := a[ j1     ] + a[ j3 + 1 ];
        x3i := a[ j1 + 1 ] - a[ j3     ];

        y0r := wd1i * x0r - wd1r * x0i;
        y0i := wd1i * x0i + wd1r * x0r;
        y2r := wk1i * x2r - wk1r * x2i;
        y2i := wk1i * x2i + wk1r * x2r;
        a[ j0     ] := y0r + y2r;
        a[ j0 + 1 ] := y0i + y2i;
        a[ j1     ] := y0r - y2r;
        a[ j1 + 1 ] := y0i - y2i;

        y0r := wd3i * x1r + wd3r * x1i;
        y0i := wd3i * x1i - wd3r * x1r;
        y2r := wk3i * x3r + wk3r * x3i;
        y2i := wk3i * x3i - wk3r * x3r;
        a[ j2     ] := y0r + y2r;
        a[ j2 + 1 ] := y0i + y2i;
        a[ j3     ] := y0r - y2r;
        a[ j3 + 1 ] := y0i - y2i;

        Inc( j, 2 )
    end;

    wk1r := w[ m ];
    wk1i := w[ m + 1 ];

    j0 := mh;
    j1 := j0 + m;
    j2 := j1 + m;
    j3 := j2 + m;

    x0r := a[ j0     ] - a[ j2 + 1 ];
    x0i := a[ j0 + 1 ] + a[ j2     ];
    x1r := a[ j0     ] + a[ j2 + 1 ];
    x1i := a[ j0 + 1 ] - a[ j2     ];
    x2r := a[ j1     ] - a[ j3 + 1 ];
    x2i := a[ j1 + 1 ] + a[ j3     ];
    x3r := a[ j1     ] + a[ j3 + 1 ];
    x3i := a[ j1 + 1 ] - a[ j3     ];
    y0r := wk1r * x0r - wk1i * x0i;
    y0i := wk1r * x0i + wk1i * x0r;
    y2r := wk1i * x2r - wk1r * x2i;
    y2i := wk1i * x2i + wk1r * x2r;
    a[ j0     ] := y0r + y2r;
    a[ j0 + 1 ] := y0i + y2i;
    a[ j1     ] := y0r - y2r;
    a[ j1 + 1 ] := y0i - y2i;

    y0r := wk1i * x1r - wk1r * x1i; y0i := wk1i * x1i + wk1r * x1r;
    y2r := wk1r * x3r - wk1i * x3i; y2i := wk1r * x3i + wk1i * x3r;
    a[ j2     ] := y0r - y2r;
    a[ j2 + 1 ] := y0i - y2i;
    a[ j3     ] := y0r + y2r;
    a[ j3 + 1 ] := y0i + y2i;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

function TDiscreteTrans.cfttree( const n,j,k:Integer ) :Integer;
var
   i, isplt, m :Integer;
begin
     if k and 3 <> 0 then
     begin
          isplt := k and 1;

          if isplt <> 0 then cftmdl1( n, _TEMP[ j - n ], _W[ _NW - ( n shr 1 ) ] )
                        else cftmdl2( n, _TEMP[ j - n ], _W[ _NW -   n         ] )
     end
     else
     begin
          m := n;
          i := k;
          while i and 3 = 0 do
          begin
               m := m shl 2;
               i := i shr 2;
          end;

          isplt := i and 1;

          if isplt <> 0 then
          begin
               while m > 128 do
               begin
                    cftmdl1( m, _TEMP[ j - m ], _W[ _NW - ( m shr 1 ) ]);

                    m := m shr 2
               end
          end
          else
          begin
               while m > 128 do
               begin
                    cftmdl2( m, _TEMP[ j - m ], _W[ _NW - m ] );

                    m := m shr 2
               end
          end
     end;

     Result := isplt
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftf161( var a:array of Double; var w:array of Double );
var
   wn4r, wk1r, wk1i,
   x00r, x00i, x01r, x01i, x02r, x02i, x03r, x03i,
   y00r, y00i, y01r, y01i, y02r, y02i, y03r, y03i,
   y04r, y04i, y05r, y05i, y06r, y06i, y07r, y07i,
   y08r, y08i, y09r, y09i, y10r, y10i, y11r, y11i,
   y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i :Double;
begin
     wn4r := w[ 1 ];
     wk1r := w[ 2 ];
     wk1i := w[ 3 ];

     x00r := a[ 00 ] + a[ 16 ]; x00i := a[ 01 ] + a[ 17 ];
     x01r := a[ 00 ] - a[ 16 ]; x01i := a[ 01 ] - a[ 17 ];
     x02r := a[ 08 ] + a[ 24 ]; x02i := a[ 09 ] + a[ 25 ];
     x03r := a[ 08 ] - a[ 24 ]; x03i := a[ 09 ] - a[ 25 ];
     y00r := x00r + x02r;       y00i := x00i + x02i;
     y04r := x00r - x02r;       y04i := x00i - x02i;
     y08r := x01r - x03i;       y08i := x01i + x03r;
     y12r := x01r + x03i;       y12i := x01i - x03r;
     x00r := a[ 02 ] + a[ 18 ];
     x00i := a[ 03 ] + a[ 19 ];
     x01r := a[ 02 ] - a[ 18 ];
     x01i := a[ 03 ] - a[ 19 ];
     x02r := a[ 10 ] + a[ 26 ];
     x02i := a[ 11 ] + a[ 27 ];
     x03r := a[ 10 ] - a[ 26 ];
     x03i := a[ 11 ] - a[ 27 ];

     y01r := x00r + x02r;               y01i := x00i + x02i;
     y05r := x00r - x02r;               y05i := x00i - x02i;
     x00r := x01r - x03i;               x00i := x01i + x03r;
     y09r := wk1r * x00r - wk1i * x00i; y09i := wk1r * x00i + wk1i * x00r;
     x00r := x01r + x03i;               x00i := x01i - x03r;
     y13r := wk1i * x00r - wk1r * x00i; y13i := wk1i * x00i + wk1r * x00r;
     x00r := a[ 4 ] + a[ 20 ];
     x00i := a[ 5 ] + a[ 21 ];
     x01r := a[ 4 ] - a[ 20 ];
     x01i := a[ 5 ] - a[ 21 ];
     x02r := a[ 12 ] + a[ 28 ];
     x02i := a[ 13 ] + a[ 29 ];
     x03r := a[ 12 ] - a[ 28 ];
     x03i := a[ 13 ] - a[ 29 ];

     y02r := x00r + x02r;              y02i := x00i + x02i;
     y06r := x00r - x02r;              y06i := x00i - x02i;
     x00r := x01r - x03i;              x00i := x01i + x03r;
     y10r := wn4r * (  x00r - x00i  ); y10i := wn4r * (  x00i + x00r  );
     x00r := x01r + x03i;              x00i := x01i - x03r;
     y14r := wn4r * (  x00r + x00i  ); y14i := wn4r * (  x00i - x00r  );
     x00r := a[ 06 ] + a[ 22 ];
     x00i := a[ 07 ] + a[ 23 ];
     x01r := a[ 06 ] - a[ 22 ];
     x01i := a[ 07 ] - a[ 23 ];
     x02r := a[ 14 ] + a[ 30 ];
     x02i := a[ 15 ] + a[ 31 ];
     x03r := a[ 14 ] - a[ 30 ];
     x03i := a[ 15 ] - a[ 31 ];

     y03r := x00r + x02r;               y03i := x00i + x02i;
     y07r := x00r - x02r;               y07i := x00i - x02i;
     x00r := x01r - x03i;               x00i := x01i + x03r;
     y11r := wk1i * x00r - wk1r * x00i; y11i := wk1i * x00i + wk1r * x00r;
     x00r := x01r + x03i;               x00i := x01i - x03r;
     y15r := wk1r * x00r - wk1i * x00i; y15i := wk1r * x00i + wk1i * x00r;
     x00r := y12r - y14r;               x00i := y12i - y14i;
     x01r := y12r + y14r;               x01i := y12i + y14i;
     x02r := y13r - y15r;               x02i := y13i - y15i;
     x03r := y13r + y15r;               x03i := y13i + y15i;
     a[ 24 ] := x00r + x02r;
     a[ 25 ] := x00i + x02i;
     a[ 26 ] := x00r - x02r;
     a[ 27 ] := x00i - x02i;
     a[ 28 ] := x01r - x03i;
     a[ 29 ] := x01i + x03r;
     a[ 30 ] := x01r + x03i;
     a[ 31 ] := x01i - x03r;

     x00r := y08r + y10r; x00i := y08i + y10i;
     x01r := y08r - y10r; x01i := y08i - y10i;
     x02r := y09r + y11r; x02i := y09i + y11i;
     x03r := y09r - y11r; x03i := y09i - y11i;
     a[ 16 ] := x00r + x02r;
     a[ 17 ] := x00i + x02i;
     a[ 18 ] := x00r - x02r;
     a[ 19 ] := x00i - x02i;
     a[ 20 ] := x01r - x03i;
     a[ 21 ] := x01i + x03r;
     a[ 22 ] := x01r + x03i;
     a[ 23 ] := x01i - x03r;

     x00r := y05r - y07i;            x00i := y05i + y07r;
     x02r := wn4r * ( x00r - x00i ); x02i := wn4r * ( x00i + x00r );
     x00r := y05r + y07i;            x00i := y05i - y07r;
     x03r := wn4r * ( x00r - x00i ); x03i := wn4r * ( x00i + x00r );
     x00r := y04r - y06i;            x00i := y04i + y06r;
     x01r := y04r + y06i;            x01i := y04i - y06r;
     a[ 08 ] := x00r + x02r;
     a[ 09 ] := x00i + x02i;
     a[ 10 ] := x00r - x02r;
     a[ 11 ] := x00i - x02i;
     a[ 12 ] := x01r - x03i;
     a[ 13 ] := x01i + x03r;
     a[ 14 ] := x01r + x03i;
     a[ 15 ] := x01i - x03r;

     x00r := y00r + y02r; x00i := y00i + y02i;
     x01r := y00r - y02r; x01i := y00i - y02i;
     x02r := y01r + y03r; x02i := y01i + y03i;
     x03r := y01r - y03r; x03i := y01i - y03i;
     a[ 00 ] := x00r + x02r;
     a[ 01 ] := x00i + x02i;
     a[ 02 ] := x00r - x02r;
     a[ 03 ] := x00i - x02i;
     a[ 04 ] := x01r - x03i;
     a[ 05 ] := x01i + x03r;
     a[ 06 ] := x01r + x03i;
     a[ 07 ] := x01i - x03r;
end;

procedure TDiscreteTrans.cftf162( var a:array of Double; var w:array of Double );
var
   wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i,
   x0r, x0i, x1r, x1i, x2r, x2i,
   y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
   y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i,
   y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i,
   y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i :Double;
begin
     wn4r := +w[ 1 ];
     wk1r := +w[ 4 ];
     wk1i := +w[ 5 ];
     wk3r := +w[ 6 ];
     wk3i := -w[ 7 ];
     wk2r := +w[ 8 ];
     wk2i := +w[ 9 ];

    x1r := a[ 0 ] - a[ 17 ];
    x1i := a[ 1 ] + a[ 16 ];
    x0r := a[ 8 ] - a[ 25 ];
    x0i := a[ 9 ] + a[ 24 ];
    x2r := wn4r * ( x0r - x0i );
    x2i := wn4r * ( x0i + x0r );
    y0r := x1r + x2r;
    y0i := x1i + x2i;
    y4r := x1r - x2r;
    y4i := x1i - x2i;
    x1r := a[ 0 ] + a[ 17 ];
    x1i := a[ 1 ] - a[ 16 ];
    x0r := a[ 8 ] + a[ 25 ];
    x0i := a[ 9 ] - a[ 24 ];
    x2r := wn4r * ( x0r - x0i );
    x2i := wn4r * ( x0i + x0r );
    y8r := x1r - x2i;
    y8i := x1i + x2r;
    y12r := x1r + x2i;
    y12i := x1i - x2r;
    x0r := a[ 2 ] - a[ 19 ];
    x0i := a[ 3 ] + a[ 18 ];
    x1r := wk1r * x0r - wk1i * x0i;
    x1i := wk1r * x0i + wk1i * x0r;
    x0r := a[ 10 ] - a[ 27 ];
    x0i := a[ 11 ] + a[ 26 ];
    x2r := wk3i * x0r - wk3r * x0i;
    x2i := wk3i * x0i + wk3r * x0r;
    y1r := x1r + x2r;
    y1i := x1i + x2i;
    y5r := x1r - x2r;
    y5i := x1i - x2i;
    x0r := a[ 2 ] + a[ 19 ];
    x0i := a[ 3 ] - a[ 18 ];
    x1r := wk3r * x0r - wk3i * x0i;
    x1i := wk3r * x0i + wk3i * x0r;
    x0r := a[ 10 ] + a[ 27 ];
    x0i := a[ 11 ] - a[ 26 ];
    x2r := wk1r * x0r + wk1i * x0i;
    x2i := wk1r * x0i - wk1i * x0r;
    y9r := x1r - x2r;
    y9i := x1i - x2i;
    y13r := x1r + x2r;
    y13i := x1i + x2i;
    x0r := a[ 4 ] - a[ 21 ];
    x0i := a[ 5 ] + a[ 20 ];
    x1r := wk2r * x0r - wk2i * x0i;
    x1i := wk2r * x0i + wk2i * x0r;
    x0r := a[ 12 ] - a[ 29 ];
    x0i := a[ 13 ] + a[ 28 ];
    x2r := wk2i * x0r - wk2r * x0i;
    x2i := wk2i * x0i + wk2r * x0r;
    y2r := x1r + x2r;
    y2i := x1i + x2i;
    y6r := x1r - x2r;
    y6i := x1i - x2i;
    x0r := a[ 4 ] + a[ 21 ];
    x0i := a[ 5 ] - a[ 20 ];
    x1r := wk2i * x0r - wk2r * x0i;
    x1i := wk2i * x0i + wk2r * x0r;
    x0r := a[ 12 ] + a[ 29 ];
    x0i := a[ 13 ] - a[ 28 ];
    x2r := wk2r * x0r - wk2i * x0i;
    x2i := wk2r * x0i + wk2i * x0r;
    y10r := x1r - x2r;
    y10i := x1i - x2i;
    y14r := x1r + x2r;
    y14i := x1i + x2i;
    x0r := a[ 6 ] - a[ 23 ];
    x0i := a[ 7 ] + a[ 22 ];
    x1r := wk3r * x0r - wk3i * x0i;
    x1i := wk3r * x0i + wk3i * x0r;
    x0r := a[ 14 ] - a[ 31 ];
    x0i := a[ 15 ] + a[ 30 ];
    x2r := wk1i * x0r - wk1r * x0i;
    x2i := wk1i * x0i + wk1r * x0r;
    y3r := x1r + x2r;
    y3i := x1i + x2i;
    y7r := x1r - x2r;
    y7i := x1i - x2i;
    x0r := a[ 6 ] + a[ 23 ];
    x0i := a[ 7 ] - a[ 22 ];
    x1r := wk1i * x0r + wk1r * x0i;
    x1i := wk1i * x0i - wk1r * x0r;
    x0r := a[ 14 ] + a[ 31 ];
    x0i := a[ 15 ] - a[ 30 ];
    x2r := wk3i * x0r - wk3r * x0i;
    x2i := wk3i * x0i + wk3r * x0r;
    y11r := x1r + x2r;
    y11i := x1i + x2i;
    y15r := x1r - x2r;
    y15i := x1i - x2i;
    x1r := y0r + y2r;
    x1i := y0i + y2i;
    x2r := y1r + y3r;
    x2i := y1i + y3i;
    a[ 0 ] := x1r + x2r;
    a[ 1 ] := x1i + x2i;
    a[ 2 ] := x1r - x2r;
    a[ 3 ] := x1i - x2i;
    x1r := y0r - y2r;
    x1i := y0i - y2i;
    x2r := y1r - y3r;
    x2i := y1i - y3i;
    a[ 4 ] := x1r - x2i;
    a[ 5 ] := x1i + x2r;
    a[ 6 ] := x1r + x2i;
    a[ 7 ] := x1i - x2r;
    x1r := y4r - y6i;
    x1i := y4i + y6r;
    x0r := y5r - y7i;
    x0i := y5i + y7r;
    x2r := wn4r * ( x0r - x0i );
    x2i := wn4r * ( x0i + x0r );
    a[ 8 ] := x1r + x2r;
    a[ 9 ] := x1i + x2i;
    a[ 10 ] := x1r - x2r;
    a[ 11 ] := x1i - x2i;
    x1r := y4r + y6i;
    x1i := y4i - y6r;
    x0r := y5r + y7i;
    x0i := y5i - y7r;
    x2r := wn4r * ( x0r - x0i );
    x2i := wn4r * ( x0i + x0r );
    a[ 12 ] := x1r - x2i;
    a[ 13 ] := x1i + x2r;
    a[ 14 ] := x1r + x2i;
    a[ 15 ] := x1i - x2r;
    x1r := y8r + y10r;
    x1i := y8i + y10i;
    x2r := y9r - y11r;
    x2i := y9i - y11i;
    a[ 16 ] := x1r + x2r;
    a[ 17 ] := x1i + x2i;
    a[ 18 ] := x1r - x2r;
    a[ 19 ] := x1i - x2i;
    x1r := y8r - y10r;
    x1i := y8i - y10i;
    x2r := y9r + y11r;
    x2i := y9i + y11i;
    a[ 20 ] := x1r - x2i;
    a[ 21 ] := x1i + x2r;
    a[ 22 ] := x1r + x2i;
    a[ 23 ] := x1i - x2r;
    x1r := y12r - y14i;
    x1i := y12i + y14r;
    x0r := y13r + y15i;
    x0i := y13i - y15r;
    x2r := wn4r * ( x0r - x0i );
    x2i := wn4r * ( x0i + x0r );
    a[ 24 ] := x1r + x2r;
    a[ 25 ] := x1i + x2i;
    a[ 26 ] := x1r - x2r;
    a[ 27 ] := x1i - x2i;
    x1r := y12r + y14i;
    x1i := y12i - y14r;
    x0r := y13r - y15i;
    x0i := y13i + y15r;
    x2r := wn4r * ( x0r - x0i );
    x2i := wn4r * ( x0i + x0r );
    a[ 28 ] := x1r - x2i;
    a[ 29 ] := x1i + x2r;
    a[ 30 ] := x1r + x2i;
    a[ 31 ] := x1i - x2r;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftf081( var a:array of Double; var w:array of Double );
var
   wn4r, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i,
   y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
   y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i :Double;
begin
     wn4r := w[ 1 ];

     x0r := a[ 00 ] + a[ 08 ];    x0i := a[ 01 ] + a[ 09 ];
     x1r := a[ 00 ] - a[ 08 ];    x1i := a[ 01 ] - a[ 09 ];
     x2r := a[ 04 ] + a[ 12 ];    x2i := a[ 05 ] + a[ 13 ];
     x3r := a[ 04 ] - a[ 12 ];    x3i := a[ 05 ] - a[ 13 ];
     y0r := x0r + x2r;            y0i := x0i + x2i;
     y2r := x0r - x2r;            y2i := x0i - x2i;
     y1r := x1r - x3i;            y1i := x1i + x3r;
     y3r := x1r + x3i;            y3i := x1i - x3r;
     x0r := a[ 02 ] + a[ 10 ];    x0i := a[ 03 ] + a[ 11 ];
     x1r := a[ 02 ] - a[ 10 ];    x1i := a[ 03 ] - a[ 11 ];
     x2r := a[ 06 ] + a[ 14 ];    x2i := a[ 07 ] + a[ 15 ];
     x3r := a[ 06 ] - a[ 14 ];    x3i := a[ 07 ] - a[ 15 ];
     y4r := x0r + x2r;            y4i := x0i + x2i;
     y6r := x0r - x2r;            y6i := x0i - x2i;
     x0r := x1r - x3i;            x0i := x1i + x3r;
     x2r := x1r + x3i;            x2i := x1i - x3r;
     y5r := wn4r * ( x0r - x0i ); y5i := wn4r * ( x0r + x0i );
     y7r := wn4r * ( x2r - x2i ); y7i := wn4r * ( x2r + x2i );

     a[ 08 ] := y1r + y5r;
     a[ 09 ] := y1i + y5i;
     a[ 10 ] := y1r - y5r;
     a[ 11 ] := y1i - y5i;
     a[ 12 ] := y3r - y7i;
     a[ 13 ] := y3i + y7r;
     a[ 14 ] := y3r + y7i;
     a[ 15 ] := y3i - y7r;
     a[ 00 ] := y0r + y4r;
     a[ 01 ] := y0i + y4i;
     a[ 02 ] := y0r - y4r;
     a[ 03 ] := y0i - y4i;
     a[ 04 ] := y2r - y6i;
     a[ 05 ] := y2i + y6r;
     a[ 06 ] := y2r + y6i;
     a[ 07 ] := y2i - y6r;
end;

procedure TDiscreteTrans.cftf082( var a:array of Double; var w:array of Double );
var
   wn4r, wk1r, wk1i, x0r, x0i, x1r, x1i,
   y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i,
   y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i :Double;
begin
    wn4r := w[ 1 ];
    wk1r := w[ 2 ];
    wk1i := w[ 3 ];

    y0r := a[ 00 ] - a[ 09 ];       y0i := a[ 01 ] + a[ 08 ];
    y1r := a[ 00 ] + a[ 09 ];       y1i := a[ 01 ] - a[ 08 ];
    x0r := a[ 04 ] - a[ 13 ];       x0i := a[ 05 ] + a[ 12 ];
    y2r := wn4r * ( x0r - x0i );    y2i := wn4r * ( x0i + x0r );
    x0r := a[ 04 ] + a[ 13 ];       x0i := a[ 05 ] - a[ 12 ];
    y3r := wn4r * ( x0r - x0i );    y3i := wn4r * ( x0i + x0r );
    x0r := a[ 02 ] - a[ 11 ];       x0i := a[ 03 ] + a[ 10 ];
    y4r := wk1r * x0r - wk1i * x0i; y4i := wk1r * x0i + wk1i * x0r;
    x0r := a[ 02 ] + a[ 11 ];       x0i := a[ 03 ] - a[ 10 ];
    y5r := wk1i * x0r - wk1r * x0i; y5i := wk1i * x0i + wk1r * x0r;
    x0r := a[ 06 ] - a[ 15 ];       x0i := a[ 07 ] + a[ 14 ];
    y6r := wk1i * x0r - wk1r * x0i; y6i := wk1i * x0i + wk1r * x0r;
    x0r := a[ 06 ] + a[ 15 ];       x0i := a[ 07 ] - a[ 14 ];
    y7r := wk1r * x0r - wk1i * x0i; y7i := wk1r * x0i + wk1i * x0r;
    x0r := y0r + y2r;               x0i := y0i + y2i;
    x1r := y4r + y6r;               x1i := y4i + y6i;
    a[ 00 ] := x0r + x1r;
    a[ 01 ] := x0i + x1i;
    a[ 02 ] := x0r - x1r;
    a[ 03 ] := x0i - x1i;

    x0r := y0r - y2r; x0i := y0i - y2i;
    x1r := y4r - y6r; x1i := y4i - y6i;
    a[ 04 ] := x0r - x1i;
    a[ 05 ] := x0i + x1r;
    a[ 06 ] := x0r + x1i;
    a[ 07 ] := x0i - x1r;

    x0r := y1r - y3i; x0i := y1i + y3r;
    x1r := y5r - y7r; x1i := y5i - y7i;
    a[ 08 ] := x0r + x1r;
    a[ 09 ] := x0i + x1i;
    a[ 10 ] := x0r - x1r;
    a[ 11 ] := x0i - x1i;

    x0r := y1r + y3i; x0i := y1i - y3r;
    x1r := y5r + y7r; x1i := y5i + y7i;
    a[ 12 ] := x0r - x1i;
    a[ 13 ] := x0i + x1r;
    a[ 14 ] := x0r + x1i;
    a[ 15 ] := x0i - x1r;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftleaf( const n,isplt:Integer; var a:array of Double );
begin
     if n = 512 then
     begin
          cftmdl1( 128, a       , _W[ _NW -  64 ] );
          cftf161(      a       , _W[ _NW -   8 ] );
          cftf162(      a[  32 ], _W[ _NW -  32 ] );
          cftf161(      a[  64 ], _W[ _NW -   8 ] );
          cftf161(      a[  96 ], _W[ _NW -   8 ] );
          cftmdl2( 128, a[ 128 ], _W[ _NW - 128 ] );
          cftf161(      a[ 128 ], _W[ _NW -   8 ] );
          cftf162(      a[ 160 ], _W[ _NW -  32 ] );
          cftf161(      a[ 192 ], _W[ _NW -   8 ] );
          cftf162(      a[ 224 ], _W[ _NW -  32 ] );
          cftmdl1( 128, a[ 256 ], _W[ _NW -  64 ] );
          cftf161(      a[ 256 ], _W[ _NW -   8 ] );
          cftf162(      a[ 288 ], _W[ _NW -  32 ] );
          cftf161(      a[ 320 ], _W[ _NW -   8 ] );
          cftf161(      a[ 352 ], _W[ _NW -   8 ] );

          if isplt <> 0 then
          begin
               cftmdl1( 128, a[ 384 ], _W[ _NW - 64 ] );
               cftf161(      a[ 480 ], _W[ _NW -  8 ] );
          end
          else
          begin
               cftmdl2( 128, a[ 384 ], _W[ _NW - 128 ] );
               cftf162(      a[ 480 ], _W[ _NW -  32 ] );
          end;

          cftf161( a[ 384 ], _W[ _NW -  8 ] );
          cftf162( a[ 416 ], _W[ _NW - 32 ] );
          cftf161( a[ 448 ], _W[ _NW -  8 ] );
     end
     else
     begin
          cftmdl1( 64, a       , _W[ _NW - 32 ] );
          cftf081(     a       , _W[ _NW -  8 ] );
          cftf082(     a[  16 ], _W[ _NW -  8 ] );
          cftf081(     a[  32 ], _W[ _NW -  8 ] );
          cftf081(     a[  48 ], _W[ _NW -  8 ] );
          cftmdl2( 64, a[  64 ], _W[ _NW - 64 ] );
          cftf081(     a[  64 ], _W[ _NW -  8 ] );
          cftf082(     a[  80 ], _W[ _NW -  8 ] );
          cftf081(     a[  96 ], _W[ _NW -  8 ] );
          cftf082(     a[ 112 ], _W[ _NW -  8 ] );
          cftmdl1( 64, a[ 128 ], _W[ _NW - 32 ] );
          cftf081(     a[ 128 ], _W[ _NW -  8 ] );
          cftf082(     a[ 144 ], _W[ _NW -  8 ] );
          cftf081(     a[ 160 ], _W[ _NW -  8 ] );
          cftf081(     a[ 176 ], _W[ _NW -  8 ] );

          if isplt <> 0 then
          begin
               cftmdl1( 64, a[ 192 ], _W[ _NW - 32 ] );
               cftf081(     a[ 240 ], _W[ _NW -  8 ] );
          end
          else
          begin
               cftmdl2( 64, a[ 192 ], _W[ _NW - 64 ] );
               cftf082(     a[ 240 ], _W[ _NW -  8 ] );
          end;

          cftf081( a[ 192 ], _W[ _NW - 8 ] );
          cftf082( a[ 208 ], _W[ _NW - 8 ] );
          cftf081( a[ 224 ], _W[ _NW - 8 ] );
    end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftrec4;
var
   isplt, j, k, m :Integer;
begin
     m := _TempN;
     while m > 512 do
     begin
          m := m shr 2;
          cftmdl1( m, _TEMP[ _TempN - m ], _W[ _NW - ( m shr 1 ) ])
     end;

     cftleaf( m, 1, _TEMP[ _TempN - m ] );

     k := 0;
     j := _TempN - m;
     while j > 0 do
     begin
          Inc( k );

          isplt := cfttree( m, j, k );
          cftleaf( m, isplt, _TEMP[ j - m ] );

          Dec( j, m )
     end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftfx41;
begin
    if _TempN = 128 then
    begin
         cftf161( _TEMP      , _W[ _NW -  8 ] );
         cftf162( _TEMP[ 32 ], _W[ _NW - 32 ] );
         cftf161( _TEMP[ 64 ], _W[ _NW -  8 ] );
         cftf161( _TEMP[ 96 ], _W[ _NW -  8 ] );
    end
    else
    begin
         cftf081( _TEMP      , _W[ _NW -  8 ] );
         cftf082( _TEMP[ 16 ], _W[ _NW -  8 ] );
         cftf081( _TEMP[ 32 ], _W[ _NW -  8 ] );
         cftf081( _TEMP[ 48 ], _W[ _NW -  8 ] );
    end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftf040;
var
   x0r, x0i,
   x1r, x1i,
   x2r, x2i,
   x3r, x3i :Double;
begin
     x0r := _TEMP[ 0 ] + _TEMP[ 4 ];  x0i := _TEMP[ 1 ] + _TEMP[ 5 ];
     x1r := _TEMP[ 0 ] - _TEMP[ 4 ];  x1i := _TEMP[ 1 ] - _TEMP[ 5 ];
     x2r := _TEMP[ 2 ] + _TEMP[ 6 ];  x2i := _TEMP[ 3 ] + _TEMP[ 7 ];
     x3r := _TEMP[ 2 ] - _TEMP[ 6 ];  x3i := _TEMP[ 3 ] - _TEMP[ 7 ];

     _TEMP[ 0 ] := x0r + x2r;  _TEMP[ 1 ] := x0i + x2i;
     _TEMP[ 2 ] := x1r - x3i;  _TEMP[ 3 ] := x1i + x3r;
     _TEMP[ 4 ] := x0r - x2r;  _TEMP[ 5 ] := x0i - x2i;
     _TEMP[ 6 ] := x1r + x3i;  _TEMP[ 7 ] := x1i - x3r;
end;

procedure TDiscreteTrans.cftb040;
var
   x0r, x0i,
   x1r, x1i,
   x2r, x2i,
   x3r, x3i :Double;
begin
     x0r := _TEMP[ 0 ] + _TEMP[ 4 ];  x0i := _TEMP[ 1 ] + _TEMP[ 5 ];
     x1r := _TEMP[ 0 ] - _TEMP[ 4 ];  x1i := _TEMP[ 1 ] - _TEMP[ 5 ];
     x2r := _TEMP[ 2 ] + _TEMP[ 6 ];  x2i := _TEMP[ 3 ] + _TEMP[ 7 ];
     x3r := _TEMP[ 2 ] - _TEMP[ 6 ];  x3i := _TEMP[ 3 ] - _TEMP[ 7 ];

     _TEMP[ 0 ] := x0r + x2r;  _TEMP[ 1 ] := x0i + x2i;
     _TEMP[ 2 ] := x1r + x3i;  _TEMP[ 3 ] := x1i - x3r;
     _TEMP[ 4 ] := x0r - x2r;  _TEMP[ 5 ] := x0i - x2i;
     _TEMP[ 6 ] := x1r - x3i;  _TEMP[ 7 ] := x1i + x3r;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftx020;
var
   x0r, x0i :Double;
begin
            x0r := _TEMP[ 0 ] - _TEMP[ 2 ];
            x0i := _TEMP[ 1 ] - _TEMP[ 3 ];

     _TEMP[ 0 ] := _TEMP[ 0 ] + _TEMP[ 2 ];
     _TEMP[ 1 ] := _TEMP[ 1 ] + _TEMP[ 3 ];

     _TEMP[ 2 ] := x0r;
     _TEMP[ 3 ] := x0i;
end;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& protected

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX ÉAÉNÉZÉX

procedure TDiscreteTrans.SetCount( const Count_:Integer );
begin
     _Count := Count_;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.cftfsub;
begin
     if _TempN > 8 then
     begin
          if _TempN > 32 then
          begin
               cftf1st( _W[ _NW - ( _TempN shr 2  ) ] );

               if _TempN > 512 then cftrec4
                               else
               if _TempN > 128 then cftleaf( _TempN, 1, _TEMP )
                               else cftfx41;

               bitrv2;
          end
          else
          if _TempN = 32 then
          begin
               cftf161( _TEMP, _W[ _NW - 8 ] );
               bitrv216;
          end
          else
          begin
               cftf081( _TEMP, _W );
               bitrv208;
          end
     end
     else
     if _TempN = 8 then cftf040
     else
     if _TempN = 4 then cftx020
end;

procedure TDiscreteTrans.cftbsub;
begin
     if _TempN > 8 then
     begin
          if _TempN > 32 then
          begin
               cftb1st( _W[ _NW - ( _TempN shr 2  ) ] );

               if _TempN > 512 then cftrec4
                               else
               if _TempN > 128 then cftleaf( _TempN, 1, _TEMP )
                               else cftfx41;

               bitrv2conj;
          end
          else
          if _TempN = 32 then
          begin
               cftf161( _TEMP, _W[ _NW - 8 ] );

               bitrv216neg;
          end
          else
          begin
               cftf081( _TEMP, _W );

               bitrv208neg;
          end
     end
     else
     if _TempN = 8 then cftb040
     else
     if _TempN = 4 then cftx020
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.rftfsub( var c:array of Double );
var
   j, k, kk, ks, m :Integer;
   wkr, wki,
    xr,  xi,
    yr,  yi :Double;
begin
     m := _TempN shr 1;

     ks := 2 * _NC div m;
     kk := 0;
     j := 2;
     while j < m do
     begin
          k := _TempN - j;
          Inc( kk, ks );

          wkr := 0.5 - c[ _NC - kk ];
          wki :=       c[       kk ];

          xr := _TEMP[ j     ] - _TEMP[ k     ];
          xi := _TEMP[ j + 1 ] + _TEMP[ k + 1 ];

          yr := wkr * xr - wki * xi;
          yi := wkr * xi + wki * xr;

          _TEMP[ j     ] := _TEMP[ j     ] - yr;
          _TEMP[ j + 1 ] := _TEMP[ j + 1 ] - yi;
          _TEMP[ k     ] := _TEMP[ k     ] + yr;
          _TEMP[ k + 1 ] := _TEMP[ k + 1 ] - yi;

          Inc( j, 2 )
     end
end;

procedure TDiscreteTrans.rftbsub( var c:array of Double );
var
   j, k, kk, ks, m :Integer;
   wkr, wki,
    xr,  xi,
    yr,  yi :Double;
begin
     m := _TempN shr 1;

     ks := 2 * _NC div m;
     kk := 0;
     j := 2;
     while j < m do
     begin
          k  := _TempN - j;
          kk := kk + ks;

          wkr := 0.5 - c[ _NC - kk ];
          wki :=       c[       kk ];

          xr := _TEMP[ j     ] - _TEMP[ k     ];
          xi := _TEMP[ j + 1 ] + _TEMP[ k + 1 ];

          yr := wkr * xr + wki * xi;
          yi := wkr * xi - wki * xr;

          _TEMP[ j     ] := _TEMP[ j     ] - yr;
          _TEMP[ j + 1 ] := _TEMP[ j + 1 ] - yi;
          _TEMP[ k     ] := _TEMP[ k     ] + yr;
          _TEMP[ k + 1 ] := _TEMP[ k + 1 ] - yi;

          Inc( j, 2 )
     end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.Normalize;
var
   N :Integer;
   P :PDouble;
begin
     P := @_TEMP[ 0 ];
     for N := 1 to _TempN do
     begin
          P^ := _NormW * P^; Inc( P )
     end
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TDiscreteTrans.MakeTableW;
begin
     _NW := _IP[ 0 ];

     if _TempN > _NW shl 2 then
     begin
          _NW := _TempN shr 2;

          makewt
     end
end;

procedure TDiscreteTrans.MakeTableC;
begin
     _NC := _IP[ 1 ];

     if _TempN > _NC shl 2 then
     begin
          _NC := _TempN shr 2;

          makect( _W[ _NW ] )
     end
end;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& public

constructor TDiscreteTrans.Create;
begin
     inherited;

end;

destructor TDiscreteTrans.Destroy;
begin

     inherited;
end;


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ÅyÉãÅ[É`ÉìÅz

//################################################################################################## Å†

initialization //############################################################################ èâä˙âª

finalization //############################################################################## èIóπâª

end. //############################################################################################# Å°
