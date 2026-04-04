---
active: true
iteration: 2
session_id: 
max_iterations: 0
completion_promise: null
started_at: "2026-04-04T00:56:14Z"
---

ssh gb10 에서 ~/data/genomes 에서 필요한 벼 레퍼런스 데이터를 다운 받는다. SRA를 다운받고 필요한 파일을 다운 받기 위해 micromamba를 이용하여 comprehensive한 환경을 만든다. 환경이 완료되면 readme.md claude.md 에 적어둔다.  만들어진 python 파이프라인을 통해서 벼 SRA 를 테스트 한다. 테스트 후에는 visualization 뿐 아니라 위치 및 필요한 모든 내용을 results.md 에 적어둔다. 성공적으로 끝났다면 x토마톨을 이용하다.   https://drive.google.com/drive/folders/1_sfaeFYIQ8bikSO9aJMDgR83DjYZaJ3r 토마토는 micorotom ep데이터를 이용한다. https://pmc.ncbi.nlm.nih.gov/articles/PMC11897730/ 토마토도 런을 이용해서 위치를 찾는다 둘다 성공적일경우 15x 10x 5x 3x 등을 reference genome size에 의해 계산해서 확인한다. 성공적이지 못할경우 코드 디버깅을 통해서 성공적으로 만든다. 매번 다시 돌리지 않고 있는 데이터를 최대한 이용한다. 마지막에 성공일 경우 리팩토링해서 이전 결과를 지우고 다시 돌린다. 작업에 문제가 있으면 디버깅해서 정상적으로 될 때까지 계속 수정해 완료시 done을 표시한다.  --max-iteration 350 --complete-promise done
