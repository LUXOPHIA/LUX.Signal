unit LUX.DiscreteTrans.Fourier.D1;

interface //######################################################################################## ■

uses LUX, LUX.D1, LUX.DiscreteTrans.D1, LUX.Complex;

type //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【型】

     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【レコード】

     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【クラス】

     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TComplexDFT

     TComplexDFT = class( TDiscreteTrans )
     private
     protected
       ///// フィールド
       _Wave :array of TDoubleC;
       _Freq :array of TDoubleC;
       ///// アクセス
       procedure SetCount( const Count_:Integer ); override;
       function GetWave( const I_:Integer ) :TDoubleC;
       procedure SetWave( const I_:Integer; const W_:TDoubleC );
       function GetFreq( const I_:Integer ) :TDoubleC;
       procedure SetFreq( const I_:Integer; const F_:TDoubleC );
     public
       constructor Create;
       destructor Destroy; override;
       ///// プロパティ
       property Wave [ const I_:Integer ] :TDoubleC read GetWave write SetWave;
       property Freq [ const I_:Integer ] :TDoubleC read GetFreq write SetFreq;
       ///// メソッド
       procedure TransWF; override;
       procedure TransFW; override;
     end;

     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRealDFT

     TRealDFT = class( TDiscreteTrans )
     private
     protected
       _WaveIter :TIter1D<Double>;
       _FreqIter :TIter1D<TDoubleC>;
       ///// アクセス
       procedure SetCount( const Count_:Integer ); override;
     public
       constructor Create;
       destructor Destroy; override;
       ///// プロパティ
       property WaveIter :TIter1D<Double>   read _WaveIter write _WaveIter;
       property FreqIter :TIter1D<TDoubleC> read _FreqIter write _FreqIter;
       ///// メソッド
       procedure TransWF; override;
       procedure TransFW; override;
     end;

//const //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【定数】

//var //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【変数】

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【ルーチン】

implementation //################################################################################### ■

uses Math;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【レコード】

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【クラス】

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TComplexDFT

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& private

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& protected

procedure TComplexDFT.SetCount( const Count_:Integer );
begin
     inherited;

     _TempN := 2 * _Count;

     SetLength( _TEMP, _TempN                                        );
     SetLength( _IP  , 2 + 1 shl Floor( log2( _Count + 0.5 ) ) div 2 );
     SetLength( _W   , _Count div 2                                  );

     _NormW := Roo2( 1 / _Count );

     MakeTableW;

     SetLength( _Wave, _Count );
     SetLength( _Freq, _Count );
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

function TComplexDFT.GetWave( const I_:Integer ) :TDoubleC;
begin
     Result := _Wave[ I_ ]
end;

procedure TComplexDFT.SetWave( const I_:Integer; const W_:TDoubleC );
begin
     _Wave[ I_ ] := W_;
end;

function TComplexDFT.GetFreq( const I_:Integer ) :TDoubleC;
begin
     Result := _Freq[ I_ ]
end;

procedure TComplexDFT.SetFreq( const I_:Integer; const F_:TDoubleC );
begin
     _Freq[ I_ ] := F_;
end;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& public

constructor TComplexDFT.Create;
begin
     inherited;

     Count := 2;
end;

destructor TComplexDFT.Destroy;
begin

     inherited;
end;

////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TComplexDFT.TransWF;
begin
     System.Move( _Wave[ 0 ], _TEMP[ 0 ], SizeOf( TDoubleC ) * _Count );

     cftfsub;

     Normalize;

     System.Move( _TEMP[ 0 ], _Freq[ 0 ], SizeOf( TDoubleC ) * _Count );
end;

procedure TComplexDFT.TransFW;
begin
     System.Move( _Freq[ 0 ], _TEMP[ 0 ], SizeOf( TDoubleC ) * _Count );

     Normalize;

     cftbsub;

     System.Move( _TEMP[ 0 ], _Wave[ 0 ], SizeOf( TDoubleC ) * _Count );
end;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRealDFT

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& private

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& protected

procedure TRealDFT.SetCount( const Count_:Integer );
begin
     inherited;

     _TempN := _Count;

     SetLength( _TEMP, _TempN                                            );
     SetLength( _IP  , 2 + 1 shl Floor( log2( _Count / 2 + 0.5 ) ) div 2 );
     SetLength( _W   , _Count div 2                                      );

     _NormW := 2 / _Count;

     MakeTableW;
     MakeTableC;
end;
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& public

constructor TRealDFT.Create;
begin
     inherited;

     Count := 2;
end;

destructor TRealDFT.Destroy;
begin

     inherited;
end;


////////////////////////////////////////////////////////////////////////////////////////////////////

procedure TRealDFT.TransWF;
var
   X :Integer;
   xi :Double;
   C :TDoubleC;
begin
     with _WaveIter do
     begin
          GoHead;

          for X := 0 to Self._Count - 1 do
          begin
               _TEMP[ X ] := Value;  GoNext
          end
     end;

     if _Count > 4 then
     begin
          cftfsub;

          rftfsub( _W[ _NW ] );
     end
     else
     if _Count = 4 then cftfsub;

     xi         := _TEMP[ 0 ] - _TEMP[ 1 ];
     _TEMP[ 0 ] := _TEMP[ 0 ] + _TEMP[ 1 ];
     _TEMP[ 1 ] := xi;

     _TEMP[ _Count - 1 ] := _TEMP[ 1 ];
     _TEMP[  1         ] := 0;

     Normalize;

     with _FreqIter do
     begin
          GoHead;

          for X := 0 to Self._Count div 2 - 1 do
          begin
               with C do
               begin
                    R := _TEMP[ 2 * X + 0 ];
                    I := _TEMP[ 2 * X + 1 ];
               end;

               Value := C;  GoNext
          end
     end
end;

procedure TRealDFT.TransFW;
var
   X :Integer;
begin
     with _FreqIter do
     begin
          GoHead;

          for X := 0 to Self._Count div 2 - 1 do
          begin
               with Value do
               begin
                    _TEMP[ 2 * X + 0 ] := R;
                    _TEMP[ 2 * X + 1 ] := I;
               end;

               GoNext
          end
     end;

     //Normalize;

     _TEMP[  1         ] := _TEMP[ _Count - 1 ];
     _TEMP[ _Count - 1 ] := 0;

     _TEMP[ 1 ] := 0.5 * ( _TEMP[ 0 ] - _TEMP[ 1 ] );
     _TEMP[ 0 ] :=         _TEMP[ 0 ] - _TEMP[ 1 ];

     if _Count > 4 then
     begin
          rftbsub( _W[ _NW ] );

          cftbsub
     end
     else
     if _Count = 4 then cftbsub;

     with _WaveIter do
     begin
          GoHead;

          for X := 0 to Self._Count - 1 do
          begin
               Value := _TEMP[ X ];  GoNext
          end
     end
end;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【ルーチン】

//################################################################################################## □

initialization //############################################################################ 初期化

finalization //############################################################################## 終了化

end. //############################################################################################# ■
