unit LUX.Signal.Fourier.D1;

interface //#################################################################### ■

uses System.Types, System.SysUtils, System.Classes,
     System.UITypes, System.Math, System.Math.Vectors,
     FMX.Types, FMX.Graphics, FMX.Controls, FMX.Objects,
     FMX.Types3D, FMX.Controls3D, FMX.MaterialSources, FMX.Objects3D,
     LUX, LUX.D1, LUX.Complex, LUX.Complex.D1, LUX.DiscreteTrans.Fourier.D1;

type //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【型】

     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【レコード】

     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【クラス】

     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FourierTimes

     TFourierTimes = class
     private
       _RDFT  :TRealDFT; 
       _TempW :TArray<Single>;
       _TempF :TArray<TSingleC>;
       _IterW :TDoubleIter1D;
       _IterF :TDoubleCIter1D;
       ///// メソッド
       procedure InitMap;
     protected
       _MapTF    :TArray2<TSingleC>;
       _TimeN    :Integer;
       _FreqN    :Integer;
       _WindN    :Integer;
       _WaveIter :TIter1D<Single>;
       ///// アクセス
       function GetMapTF( const T_,F_:Integer ) :TSingleC;
       procedure SetTimeN( const TimeN_:Integer );
       procedure SetFreqN( const FreqN_:Integer );
       procedure SetWindN( const WindN_:Integer );
     public
       constructor Create( const TimeN_,FreqN_:Integer; const WaveIter_:TIter1D<Single> );
       destructor Destroy; override;
       ///// プロパティ
       property MapTF[ const T_,F_:Integer ] :TSingleC        read GetMapTF                     ;
       property TimeN                        :Integer         read   _TimeN    write SetTimeN   ;
       property FreqN                        :Integer         read   _FreqN    write SetFreqN   ;
       property WindN                        :Integer         read   _WindN    write SetWindN   ;
       property WaveIter                     :TIter1D<Single> read   _WaveIter write   _WaveIter;
       ///// メソッド
       procedure Analyze;
     end;

//const //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【定数】

//var //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【変数】

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【ルーチン】

implementation //############################################################### ■

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【レコード】

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【クラス】

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FourierTimes

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& private

procedure TFourierTimes.InitMap;
begin
     SetLength( _MapTF, _TimeN, _FreqN );

     SetLength( _TempW, _WindN );
     SetLength( _TempF, _FreqN );

     TArrayIter< Single >( _IterW.Iter ).Parent := _TempW;
     TArrayIter<TSingleC>( _IterF.Iter ).Parent := _TempF;

     _RDFT.Count := _WindN;
end;

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& protected

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX アクセス

function TFourierTimes.GetMapTF( const T_,F_:Integer ) :TSingleC;
begin
     Result := _MapTF[ T_, F_ ];
end;

procedure TFourierTimes.SetTimeN( const TimeN_:Integer );
begin
     _TimeN := TimeN_;

     InitMap;
end;

procedure TFourierTimes.SetFreqN( const FreqN_:Integer );
begin
     _FreqN := FreqN_;

     _WindN := 2 * _FreqN;

     InitMap;
end;

procedure TFourierTimes.SetWindN( const WindN_:Integer );
begin
     _WindN := WindN_;

     _FreqN := _WindN div 2;

     InitMap;
end;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX メソッド

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& public

constructor TFourierTimes.Create( const TimeN_,FreqN_:Integer; const WaveIter_:TIter1D<Single> );
begin
     inherited Create;

     _IterW := TDoubleIter1D .Create( TArrayIter< Single >.Create( _TempW ) );
     _IterF := TDoubleCIter1D.Create( TArrayIter<TSingleC>.Create( _TempF ) );

     _RDFT := TRealDFT.Create;
     with _RDFT do
     begin
          WaveIter := _IterW;
          FreqIter := _IterF;
     end;

     _TimeN    := TimeN_;
      FreqN    := FreqN_;
     _WaveIter := WaveIter_;
end;

destructor TFourierTimes.Destroy;
begin
     _RDFT .Free;

     _IterW.Free;
     _IterF.Free;

     inherited;
end;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX メソッド

function Blackman( const X_:Double ) :Double;
begin
     if ( 0 < X_ ) and ( X_ < 1 ) then Result := 0.42
                                               - 0.5  * Cos( Pi2 * X_ )
                                               + 0.08 * Cos( Pi4 * X_ )
                                  else Result := 0
end;

procedure TFourierTimes.Analyze;
var
   N, X, Y :Integer;
begin
     N := ( _WaveIter.Count - _WindN ) div _TimeN;

     for X := 0 to _TimeN - 1 do
     begin
          with _WaveIter do
          begin
               GoJump( N * X );
               for Y := 0 to _WindN-1 do
               begin
                    _TempW[ Y ] := Blackman( ( 0.5 + Y ) / (512-1) ) * Value;  GoNext;
               end;
          end;

          _RDFT.TransWF;

          for Y := 0 to _FreqN-1 do _MapTF[ X, Y ] := _TempF[ Y ];
     end;
end;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$【ルーチン】

//############################################################################## □

initialization //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 初期化

finalization //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 最終化

end. //######################################################################### ■
