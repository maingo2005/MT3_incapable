#include <Novice.h>

//#include "../MT3.h"
#include <imgui.h>

const char kWindowTitle[] = "MT3_03_01_Basic";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	const int kWindowWidth = 1280;
	const int kWindowHeight = 720;
	const float kAxesSize = 100.0f;
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Vector3 translates[3] = {
	  {0.2f, 1.0f, 0.0f},
	  {0.4f, 0.0f, 0.0f},
	  {0.3f, 0.0f, 0.0f},
	};
	Vector3 rotates[3] = {
	  {0.0f, 0.0f, -6.8f},
	  {0.0f, 0.0f, -1.4f},
	  {0.0f, 0.0f, 0.0f },
	};
	Vector3 scales[3] = {
	  {1.0f, 1.0f, 1.0f},
      {1.0f, 1.0f, 1.0f},
      {1.0f, 1.0f, 1.0f}
    };

	Camera camera{};
	camera.spherical.theta = 0.26f;
	camera.spherical.radius = 7.0f;

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		Update(camera);

		ImGui::Begin("Window");

		ImGui::DragFloat3("translates[0]", &translates[0].x, 0.01f);
		ImGui::DragFloat3("rotates[0]", &rotates[0].x, 0.01f);
		ImGui::DragFloat3("scales[0]", &scales[0].x, 0.01f);

		ImGui::DragFloat3("translates[1]", &translates[1].x, 0.01f);
		ImGui::DragFloat3("rotates[1]", &rotates[1].x, 0.01f);
		ImGui::DragFloat3("scales[1]", &scales[1].x, 0.01f);

		ImGui::DragFloat3("translates[2]", &translates[2].x, 0.01f);
		ImGui::DragFloat3("rotates[2]", &rotates[2].x, 0.01f);
		ImGui::DragFloat3("scales[2]", &scales[2].x, 0.01f);

		ImGui::End();

		Matrix4x4 localMatrix0 = MakeAffineMatrix(scales[0], rotates[0], translates[0]);
		Matrix4x4 localMatrix1 = MakeAffineMatrix(scales[1], rotates[1], translates[1]);
		Matrix4x4 localMatrix2 = MakeAffineMatrix(scales[2], rotates[2], translates[2]);

		Matrix4x4 worldMatrix0 = localMatrix0;
		Matrix4x4 worldMatrix1 = Multiply(localMatrix1, worldMatrix0);
		Matrix4x4 worldMatrix2 = Multiply(localMatrix2, worldMatrix1);

		Sphere sphere0{
		  .center = {worldMatrix0.m[3][0], worldMatrix0.m[3][1], worldMatrix0.m[3][2]},
		  .radius = 0.05f,
		};
		Sphere sphere1{
		  .center = {worldMatrix1.m[3][0], worldMatrix1.m[3][1], worldMatrix1.m[3][2]},
		  .radius = 0.05f,
		};
		Sphere sphere2{
		  .center = {worldMatrix2.m[3][0], worldMatrix2.m[3][1], worldMatrix2.m[3][2]},
		  .radius = 0.05f
        };

		Matrix4x4 viewMatrix = CalcViewMatrix(camera);
		Matrix4x4 projectionMatrix =
		  MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);

		// ViewportMatrixを作る
		Matrix4x4 viewportMatrix =
		  MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		Matrix4x4 viewportMatrixForAxes = MakeViewportMatrix(
		  0.0f, float(kWindowHeight) - kAxesSize, kAxesSize, kAxesSize, 0.0f, 1.0f);

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		// 描画
		DrawGrid(viewProjectionMatrix, viewportMatrix);

		DrawSphere(sphere0, viewProjectionMatrix, viewportMatrix, RED);
		DrawSphere(sphere1, viewProjectionMatrix, viewportMatrix, GREEN);
		DrawSphere(sphere2, viewProjectionMatrix, viewportMatrix, BLUE);

		DrawSegment(
		  Segment{sphere0.center, sphere1.center - sphere0.center}, viewProjectionMatrix,
		  viewportMatrix, WHITE);
		DrawSegment(
		  Segment{sphere1.center, sphere2.center - sphere1.center}, viewProjectionMatrix,
		  viewportMatrix, WHITE);

		DrawAxis(camera, viewportMatrixForAxes);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
