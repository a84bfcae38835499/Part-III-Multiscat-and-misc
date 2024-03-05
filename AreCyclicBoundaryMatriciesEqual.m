function [result] = AreCyclicBoundaryMatriciesEqual(A,B)
  if(~isequal(size(A),size(B)))
    disp(size(A))
    disp(size(B))
    error("Matricies not of equal size!")
  end
  Nsuper = length(A);
  %operations to check:
  % - rot90 (and 180,270)
  % - fliplr and flipud
  % - transpose
  %oh god wait some of these are the same. Like, a transpose is the same as
  %flipud followed by rot90. how can we generalise transformations and make
  %sure they're not repeated??
  %Hypothesis: Doing each rot, checking that, then checking its transpose,
  %then undoing the transpose and rotating again, gives the right result.
  %Then I just do that for each possible collumn/row shift!
  result = false;
  for mshift = 0:Nsuper-1
    for nshift = 0:Nsuper-1
      %Keep a constant and compare to a shifted version of B
      Bs = [B(:,2+mshift:end)  B(:,1:1+mshift)];
      Bs = [Bs(2+nshift:end,:);Bs(1:1+nshift,:)];
      Br = Bs;
      for rotN = 1:4
        Br = rot90(Br);
        if(isequal(A,Br) || isequal(A,transpose(Br)))
          result = true;
          return;
        end
        %N * N * 4 * 2
        %shift*shift*rot*transpose
      end
    end
  end
end