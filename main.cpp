#include <Novice.h>

//#include "../MT3.h"
#include <imgui.h>

const char kWindowTitle[] = "MT3_04_04_Basic";

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	struct Ball {
		Vector3 position;     // ボールの位置
		Vector3 velocity;     // ボールの速度
		Vector3 acceleration; // ボールの加速度
		float mass;           // ボールの質量
		float radius;         // ボールの半径
		unsigned int color;   // ボールの色
	};

	// ライブラリの初期化
	const int kWindowWidth = 1280;
	const int kWindowHeight = 720;
	const float kAxesSize = 100.0f;
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Camera camera{};
	camera.spherical.theta = 0.26f;
	camera.spherical.radius = 7.0f;

	Plane plane;
	plane.normal = Normalize({-0.2f, 0.9f, -0.3f});
	plane.distance = 0.0f;

	Ball ball{};
	ball.position = {0.8f, 1.2f, 0.3f};
	ball.mass = 2.0f;
	ball.radius = 0.05f;
	ball.color = WHITE;

	float e = 0.7f;
	bool start = false;
	float deltaTime = 1.0f / 60.0f;

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
		ImGui::SliderFloat3("planeNormal", &plane.normal.x, -1.0f, 1.0f);
		ImGui::SliderFloat("planeDistance", &plane.distance, -5.0f, 5.0f);
		plane.normal = Normalize(plane.normal);

		if (ImGui::Button("start")) {
			start = true;
			ball.acceleration = {0.0f, -9.8f, 0.0f};
			ball.position = {0.8f, 1.2f, 0.3f};
			ball.velocity = {0.0f, 0.0f, 0.0f};
		}

		ImGui::End();

		ball.velocity += ball.acceleration * deltaTime;
		ball.position += ball.velocity * deltaTime;

		if (IsCollision(Sphere{ball.position, ball.radius}, plane)) {
			ball.velocity = Reflect(ball.velocity, plane.normal) * e;
		}

		Matrix4x4 viewMatrix = CalcViewMatrix(camera);
		Matrix4x4 projectionMatrix =
		  MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = viewMatrix * projectionMatrix;

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

		DrawPlane(plane, viewProjectionMatrix, viewportMatrix, WHITE);

		DrawSphere(
		  Sphere{.center = ball.position, .radius = ball.radius}, viewProjectionMatrix,
		  viewportMatrix, ball.color);

		// DrawSphere(
		//   Sphere{.center = origin, .radius = 0.3f}, viewProjectionMatrix, viewportMatrix, BLUE);

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
