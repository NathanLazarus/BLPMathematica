% Validate - makes sure that numbers are not NaN or Inf

function [ bValid ] = Validate( x )

  bIsNan = any( isnan( x(:) ) ) ;
  bIsInf = any( isinf( x(:) ) ) ;
  bValid = ~bIsNan & ~bIsInf ;
  
  if ~bValid
    if bIsInf
      fprintf( 2, '---> IsInf\n' ) ;
    end
    
    if bIsNan
      fprintf( 2, '---> IsNan\n' ) ;
    end
    
% %     keyboard ;
  end